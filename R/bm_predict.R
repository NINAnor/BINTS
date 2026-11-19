# Prediction engine -----------------------------------------------------------
#
# Posterior draws of the linear predictor ("location") mu_r are assembled for
# every requested row from the monitored parameters, term by term (intercept,
# fixed effects, temporal, spatial, interaction, random effects). Everything
# downstream (expected response, predictive draws, posterior predictive
# checks) is a function of mu_r and sigma_r, using the family registry, so the
# same likelihood and truncation as in the fitted model is used throughout.
#
# Three quantities are kept distinct:
#   location   mu              location of the untruncated normal (or log-mean)
#   expected   E[y | mu, sigma]  mean of the truncated distribution
#   predictive y_rep ~ p(y | mu, sigma)   a new observation, incl. noise

bm_draw_index <- function(n_total, n_draws) {
  if (is.null(n_draws) || n_draws >= n_total) return(seq_len(n_total))
  unique(round(seq(1, n_total, length.out = n_draws)))
}

bm_label <- function(x, effect) {
  if (effect == "replication") {
    num <- suppressWarnings(as.integer(sub(".*[^0-9]([0-9]+)$", "\\1", as.character(x))))
    if (!anyNA(num)) x <- num
  }
  as.character(x)
}

bm_resolve_random_modes <- function(random, random_effects) {
  modes <- stats::setNames(rep("conditional", length(random_effects)), random_effects)
  ok <- c("conditional", "marginal", "zero")
  if (is.null(random)) return(modes)
  if (!is.character(random) || !all(random %in% ok)) {
    stop("`random` values must be from: ", paste(ok, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(names(random))) {
    if (length(random) != 1L) {
      stop("An unnamed `random` must be a single mode; use a named vector for ",
           "per-effect modes.", call. = FALSE)
    }
    modes[] <- random
  } else {
    bad <- setdiff(names(random), random_effects)
    if (length(bad)) {
      stop("`random` names not among the model's random effects (",
           paste(random_effects, collapse = ", "), "): ", paste(bad, collapse = ", "),
           ".", call. = FALSE)
    }
    modes[names(random)] <- random
  }
  modes
}

# Standardised design matrix for new data, using the transformation learned on
# the fitting data. NA covariates give NA rows.
bm_new_design <- function(fit, r, dat) {
  des <- fit$prep$design[[r]]
  if (is.null(des)) return(NULL)
  vars <- fit$spec$fixed_effects[[r]]
  miss <- setdiff(vars, names(dat))
  if (length(miss)) {
    stop("`newdata` lacks covariate column(s) needed for '", r, "': ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  old <- fit$prep$data
  xlev <- lapply(stats::setNames(vars, vars), function(v) {
    if (is.numeric(old[[v]])) NULL else levels(factor(old[[v]]))
  })
  xlev <- Filter(Negate(is.null), xlev)
  f <- stats::as.formula(paste("~", paste0("`", vars, "`", collapse = " + ")))
  mf <- tryCatch(
    stats::model.frame(f, dat, na.action = stats::na.pass, xlev = xlev),
    error = function(e) {
      stop("Covariate levels in `newdata` do not match the fitting data: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  mm <- stats::model.matrix(f, mf)[, -1, drop = FALSE]
  mm <- mm[, des$scaling$term, drop = FALSE]
  sweep(sweep(mm, 2, des$scaling$center, "-"), 2, des$scaling$scale, "/")
}

# Everything needed to evaluate the model at the requested rows
bm_pred_inputs <- function(fit, newdata, responses) {
  prep <- fit$prep
  spec <- fit$spec
  cols <- spec$columns
  has_site <- !is.null(prep$spatial)

  if (is.null(newdata)) {
    dat <- prep$data
    base <- prep$data
    rows <- prep$kept_rows
    X <- stats::setNames(lapply(responses, function(r) prep$design[[r]]$X), responses)
    year_c <- if (identical(prep$temporal$type, "linear")) prep$temporal$year_c
    season <- if (identical(prep$temporal$type, "cyclic")) prep$temporal$season
    re <- lapply(stats::setNames(spec$random_effects, spec$random_effects), function(e) {
      list(id = prep$index[[e]]$id,
           label = prep$index[[e]]$levels[prep$index[[e]]$id])
    })
    site_idx <- if (has_site) prep$spatial$site_id
    xy <- if (has_site) prep$spatial$xy[prep$spatial$site_id, , drop = FALSE]
    return(list(n = nrow(dat), dat = dat, base = base, rows = rows, X = X,
                year_c = year_c, season = season, re = re, site_idx = site_idx,
                xy = xy, from_data = TRUE))
  }

  base <- newdata
  if (inherits(newdata, "sf")) {
    geom <- sf::st_geometry(newdata)
    crs <- if (has_site) prep$spatial$options$crs else 4326
    if (is.na(sf::st_crs(geom))) sf::st_crs(geom) <- sf::st_crs(crs)
    geom <- sf::st_transform(geom, sf::st_crs(crs))
    pts <- suppressWarnings(sf::st_centroid(geom))
    cc <- sf::st_coordinates(pts)
    dat <- as.data.frame(sf::st_drop_geometry(newdata))
    if (!is.null(cols$longitude) && !cols$longitude %in% names(dat)) dat[[cols$longitude]] <- cc[, 1]
    if (!is.null(cols$latitude) && !cols$latitude %in% names(dat)) dat[[cols$latitude]] <- cc[, 2]
  } else if (is.data.frame(newdata)) {
    dat <- as.data.frame(newdata)
  } else {
    stop("`newdata` must be NULL, a data frame or an sf object.", call. = FALSE)
  }
  n <- nrow(dat)
  if (n == 0L) stop("`newdata` has no rows.", call. = FALSE)

  need <- function(col, what) {
    if (is.null(col) || !col %in% names(dat)) {
      stop("`newdata` lacks column '", col, "' needed for ", what, ".", call. = FALSE)
    }
  }
  X <- stats::setNames(lapply(responses, function(r) bm_new_design(fit, r, dat)), responses)

  year_c <- season <- NULL
  if (spec$temporal_effect == "linear") {
    need(cols$year, "the linear temporal effect")
    if (!is.numeric(dat[[cols$year]])) stop("Year column '", cols$year, "' must be numeric.", call. = FALSE)
    year_c <- dat[[cols$year]] - prep$temporal$year_center
  } else if (spec$temporal_effect == "cyclic") {
    need(cols$season, "the cyclic temporal effect")
    season <- dat[[cols$season]]
    period <- prep$temporal$cyclic$period
    if (any(season < 0 | season > period, na.rm = TRUE)) {
      stop("Seasonal values must lie within [0, ", period, "].", call. = FALSE)
    }
  }

  re <- lapply(stats::setNames(spec$random_effects, spec$random_effects), function(e) {
    col <- cols[[e]]
    if (is.null(col) || !col %in% names(dat)) {
      return(list(id = rep(NA_integer_, n), label = rep(NA_character_, n)))
    }
    lab <- bm_label(dat[[col]], e)
    list(id = match(lab, prep$index[[e]]$levels), label = lab)
  })

  site_idx <- xy <- NULL
  if (has_site) {
    lc <- cols$locality
    site_idx <- if (!is.null(lc) && lc %in% names(dat)) {
      match(as.character(dat[[lc]]), prep$spatial$site_names)
    } else rep(NA_integer_, n)
    xcol <- cols$longitude
    ycol <- cols$latitude
    if (all(c(xcol, ycol) %in% names(dat))) {
      xy <- cbind(dat[[xcol]], dat[[ycol]])
    } else {
      if (anyNA(site_idx)) {
        stop("`newdata` needs coordinate columns '", xcol, "' and '", ycol,
             "' (or an sf geometry) for locations not in the fitting data.",
             call. = FALSE)
      }
      xy <- prep$spatial$xy[site_idx, , drop = FALSE]
    }
    # Known sites use their fitted coordinates
    known <- !is.na(site_idx)
    xy[known, ] <- prep$spatial$xy[site_idx[known], , drop = FALSE]
  }
  list(n = n, dat = dat, base = base, rows = seq_len(n), X = X, year_c = year_c,
       season = season, re = re, site_idx = site_idx, xy = xy, from_data = FALSE)
}

# Draws of a spatial field (tag "sp" or "spt") of response r at each row
bm_field_at_rows <- function(fit, r, tag, inp, mat) {
  sp <- fit$prep$spatial
  D <- nrow(mat)
  out <- matrix(NA_real_, D, inp$n)
  known <- !is.na(inp$site_idx)
  if (any(known)) {
    if (sp$type == "autocorrelation") {
      out[, known] <- mat[, paste0("w_", tag, "_", r, "[", inp$site_idx[known], "]"), drop = FALSE]
    } else {
      us <- unique(inp$site_idx[known])
      fld <- bm_spatial_draws(fit, r, sp$xy[us, , drop = FALSE], tag = tag, mat = mat)
      out[, known] <- fld[, match(inp$site_idx[known], us), drop = FALSE]
    }
  }
  unk <- !known & stats::complete.cases(inp$xy)
  if (any(unk)) {
    key <- paste(inp$xy[unk, 1], inp$xy[unk, 2])
    uk <- !duplicated(key)
    fld <- bm_spatial_draws(fit, r, inp$xy[unk, , drop = FALSE][uk, , drop = FALSE],
                            simulate = TRUE, tag = tag, mat = mat)
    out[, unk] <- fld[, match(key, key[uk]), drop = FALSE]
  }
  out
}

# Draws of new random-effect values for G unseen groups, respecting how the
# effect is shared between responses: independent, a shared latent process
# (second = lambda * reference), or correlated.
bm_sim_new_effects <- function(fit, e, G, mat) {
  spec <- fit$spec
  resp <- spec$responses
  sh <- spec$sharing
  D <- nrow(mat)
  sdr <- function(r) mat[, paste0("sd_", e, "_", r)]
  Z <- function() matrix(stats::rnorm(D * G), D, G)
  out <- list()
  Z1 <- Z()
  out[[resp[1]]] <- sdr(resp[1]) * Z1
  if (length(resp) > 1L) {
    r2 <- resp[2]
    if (isTRUE(sh$active) && e %in% sh$random) {
      lam <- if (isTRUE(sh$loadings)) mat[, paste0("lambda_re_", e, "_", r2)] else 1
      out[[r2]] <- lam * out[[resp[1]]]
    } else if (isTRUE(sh$active) && e %in% sh$correlated) {
      rho <- mat[, paste0("rho_", e, "_", r2)]
      out[[r2]] <- sdr(r2) * (rho * Z1 + sqrt(1 - rho^2) * Z())
    } else {
      out[[r2]] <- sdr(r2) * Z()
    }
  }
  out
}

# Contribution of random effect `e` to mu, for every response
bm_re_contrib <- function(fit, e, inp, mode, new_levels, mat, responses) {
  if (mode == "zero") return(NULL)
  D <- nrow(mat)
  n <- inp$n
  info <- inp$re[[e]]
  nlev <- length(fit$prep$index[[e]]$levels)
  known <- !is.na(info$id)
  sim_rows <- if (mode == "marginal") rep(TRUE, n) else !known
  if (mode == "conditional" && any(!known) && new_levels == "error") {
    lab <- unique(info$label[!known])
    stop("`newdata` has ", e, " level(s) not in the fitting data (",
         paste(utils::head(lab, 5), collapse = ", "),
         "); use `new_levels = \"simulate\"` or \"zero\".", call. = FALSE)
  }
  out <- stats::setNames(vector("list", length(responses)), responses)
  for (r in responses) out[[r]] <- matrix(0, D, n)
  if (mode == "conditional" && any(known)) {
    for (r in responses) {
      U <- mat[, paste0("u_", e, "_", r, "[", seq_len(nlev), "]"), drop = FALSE]
      out[[r]][, known] <- U[, info$id[known], drop = FALSE]
    }
  }
  if (any(sim_rows) && !(mode == "conditional" && new_levels == "zero")) {
    # Rows sharing a group label share a draw; unlabelled rows are independent
    lab <- info$label
    key <- ifelse(is.na(lab), paste0("..row", seq_len(n)), lab)
    key <- key[sim_rows]
    uk <- unique(key)
    sims <- bm_sim_new_effects(fit, e, length(uk), mat)
    gi <- match(key, uk)
    for (r in responses) out[[r]][, sim_rows] <- sims[[r]][, gi, drop = FALSE]
  }
  out
}

# Posterior draws of mu and sigma at the requested rows
bm_predict_core <- function(fit, newdata, responses, random, new_levels, spatial,
                            n_draws, seed) {
  spec <- fit$spec
  if (is.null(responses)) responses <- spec$responses
  if (!is.character(responses) || !all(responses %in% spec$responses)) {
    stop("`responses` must be among the fitted responses: ",
         paste(spec$responses, collapse = ", "), ".", call. = FALSE)
  }
  modes <- bm_resolve_random_modes(random, spec$random_effects)
  inp <- bm_pred_inputs(fit, newdata, responses)

  all_mat <- as.matrix(fit$samples)
  if (!is.null(n_draws) && (!is.numeric(n_draws) || length(n_draws) != 1L || n_draws < 1)) {
    stop("`n_draws` must be a positive integer or NULL.", call. = FALSE)
  }
  idx <- bm_draw_index(nrow(all_mat), n_draws)
  mat <- all_mat[idx, , drop = FALSE]
  D <- nrow(mat)

  run <- function() {
    mu <- stats::setNames(vector("list", length(responses)), responses)
    for (r in responses) {
      m <- matrix(mat[, paste0("beta0_", r)], D, inp$n)
      X <- inp$X[[r]]
      if (!is.null(X)) {
        B <- mat[, paste0("beta_", r, "[", seq_len(ncol(X)), "]"), drop = FALSE]
        m <- m + B %*% t(X)
      }
      fc <- NULL
      if (spec$temporal_effect == "linear") {
        m <- m + outer(mat[, paste0("beta_trend_", r)], inp$year_c)
      } else if (spec$temporal_effect == "cyclic") {
        ok <- !is.na(inp$season)
        fc <- matrix(NA_real_, D, inp$n)
        fc[, ok] <- bm_cyclic_draws(fit, r, inp$season[ok], mat = mat)
        m <- m + fc
      }
      if (spec$spatial_effect != "none" && spatial == "conditional") {
        m <- m + bm_field_at_rows(fit, r, "sp", inp, mat)
        if (isTRUE(spec$spatiotemporal_interaction)) {
          g <- bm_field_at_rows(fit, r, "spt", inp, mat)
          m <- m + if (spec$temporal_effect == "linear") {
            sweep(g, 2, inp$year_c, "*")
          } else g * fc
        }
      }
      mu[[r]] <- m
    }
    for (e in spec$random_effects) {
      con <- bm_re_contrib(fit, e, inp, modes[[e]], new_levels, mat, responses)
      if (!is.null(con)) for (r in responses) mu[[r]] <- mu[[r]] + con[[r]]
    }
    mu
  }
  mu <- if (is.null(seed)) run() else bm_with_seed(seed, run())
  # second parameter of each family: residual sd `sigma`, or size `phi`
  sigma <- stats::setNames(lapply(responses, function(r) {
    p2 <- bm_families()[[spec$family[[r]]]]$par2
    if (is.null(p2)) rep(NA_real_, nrow(mat)) else mat[, paste0(p2, "_", r)]
  }), responses)
  list(mu = mu, sigma = sigma, inputs = inp, draw_index = idx, responses = responses)
}

bm_summarise_draws <- function(M, level) {
  a <- (1 - level) / 2
  n <- ncol(M)
  out <- data.frame(mean = rep(NA_real_, n), median = NA_real_, sd = NA_real_,
                    lower = NA_real_, upper = NA_real_)
  ok <- colSums(is.na(M)) == 0
  if (any(ok)) {
    Mo <- M[, ok, drop = FALSE]
    q <- apply(Mo, 2, stats::quantile, probs = c(a, 0.5, 1 - a), names = FALSE)
    out$mean[ok] <- colMeans(Mo)
    out$median[ok] <- q[2, ]
    out$sd[ok] <- apply(Mo, 2, stats::sd)
    out$lower[ok] <- q[1, ]
    out$upper[ok] <- q[3, ]
  }
  out
}

#' Predict from a fitted biomass-richness model
#'
#' Posterior predictions at the original observations, at new years and
#' covariate values, at new locations or on a spatial grid, either conditional
#' on the estimated random effects or marginalised over them.
#'
#' @section What is predicted:
#' Three quantities are kept distinct, and `type` selects one:
#' * `"expected"` (default): the posterior of the **expected response**
#'   \eqn{E[y \mid \mu, \sigma]}, the mean of the truncated distribution. For the
#'   truncated normal this is
#'   \eqn{\mu + \sigma\,\phi(\mu/\sigma)/\Phi(\mu/\sigma)}, which is always
#'   positive and at least `mu` (they coincide when `mu` is many `sigma`
#'   above zero); for the log-normal it is
#'   \eqn{\exp(\mu + \sigma^2/2)}; for the Poisson and negative binomial it is
#'   \eqn{\exp(\mu)}, the mean count. Uncertainty comes from the parameters (and from
#'   any simulated random effects), not from observation noise.
#' * `"location"`: the linear predictor `mu`, i.e. the location of the
#'   *untruncated* normal (or mean of log y, or the log mean count for the
#'   Poisson and negative binomial). It can be negative and, for the truncated
#'   normal, is **not** the mean of the response.
#' * `"predictive"`: draws of a new observation \eqn{y^{rep}} from the
#'   fitted distribution (with observation noise), so the intervals are
#'   prediction intervals. These are never negative (and are integers for the
#'   count families). The two responses are
#'   conditionally independent given their linear predictors.
#'
#' @section Random effects and new levels:
#' `random` sets how random effects (and only them) enter; it is one mode for all
#' effects or a named vector such as `c(year = "marginal", locality = "zero")`:
#' * `"conditional"`: the estimated effect of each level present in the fitting
#'   data. Levels not in the fitting data (new years, new localities, new
#'   replicates, or rows where the column is absent) are handled by
#'   `new_levels`.
#' * `"marginal"`: population-level predictions that integrate over the effect:
#'   in every posterior draw each group (rows with the same label share it)
#'   gets a fresh draw from `Normal(0, sd)`; the averages over draws are
#'   Monte Carlo marginal quantities, appropriate for the nonlinear
#'   truncated mean. Rows with the column absent are independent groups.
#' * `"zero"`: the effect is set to 0, the typical group (the conditional mode
#'   of a new group, not the population average).
#'
#' `new_levels` (for `"conditional"`) is `"simulate"` (default; draw from
#' `Normal(0, sd)`, which carries the uncertainty about unseen groups),
#' `"zero"`, or `"error"`. The simulation respects the model structure: a
#' shared latent effect gives the second response `lambda` times the same
#' draw, and correlated effects are drawn jointly with correlation `rho`.
#'
#' @section Temporal and spatial terms:
#' A linear trend is extrapolated to new years; a cyclic effect needs the
#' seasonal variable within `[0, period]`. Spatial effects are evaluated at
#' known localities (matched by the locality column) from the fitted field. For
#' unsampled locations (coordinates only, or an `sf` object whose geometry is
#' used, centroids for polygons, transformed to the model's CRS), the NNGP
#' field is drawn from its conditional distribution given the `m` nearest
#' fitted sites and the posterior of its parameters; the coordinate spline is
#' evaluated at the new coordinates (deterministic, and it should not be
#' trusted outside the sampled area). `spatial = "none"` sets spatial fields
#' to zero. Draws for the two responses at new NNGP locations are made
#' independently (their marginal distributions are right, their dependence is
#' not propagated). Rows with missing values in any required input return `NA`.
#'
#' @param object A `bm_fit` object from [fit_model()].
#' @param newdata `NULL` (the observations used in the fit), a data frame, or
#'   an `sf` object, containing the columns used by the model (covariates, year or
#'   season, coordinates and/or locality; random-effect columns are optional).
#' @param responses Subset of the fitted responses (default all).
#' @param type `"expected"`, `"location"` or `"predictive"`; see above.
#' @param random,new_levels See above.
#' @param spatial `"conditional"` (default) or `"none"`.
#' @param n_draws Number of posterior draws to use (evenly spaced; `NULL` for
#'   all). Memory is `n_draws` x rows doubles per response.
#' @param level Credible (or prediction) interval level, default 0.95.
#' @param summarise If `TRUE` (default) return posterior summaries; if `FALSE`
#'   return the draws.
#' @param seed Optional seed for the simulated parts (new levels, new-location
#'   fields, predictive draws); the RNG state of the session is restored.
#' @param ... Unused.
#'
#' @return With `summarise = TRUE`, a data frame (an `sf` object if `newdata`
#'   is `sf`) in long format with one row per input row and response: the input
#'   columns, `.row` (row of the original data when `newdata = NULL`),
#'   `response`, and the posterior `mean`, `median`, `sd`, `lower` and `upper`
#'   (equal-tailed interval at `level`). The attribute `observations` records
#'   whether the rows are the observations used in the fit (`newdata = NULL`). With `summarise = FALSE`, a list of
#'   class `bm_draws` with `draws` (a list of draws x rows matrices, one per
#'   response), `data` (the input rows), `.row`, `type` and `draw_index`.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 12, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", fixed_effects = "x1",
#'                  n_chains = 2, n_iter = 2000, n_burnin = 1000, n_thin = 2,
#'                  seed = 1, parallel = FALSE)
#' head(predict(fit))
#' # population-level prediction for a new year and locality
#' nd <- data.frame(year = 2030, x1 = c(-1, 0, 1))
#' predict(fit, nd, random = "marginal", seed = 2)
#' # draws of new observations
#' str(predict(fit, nd, type = "predictive", summarise = FALSE, n_draws = 200)$draws)
#' }
#' @md
predict.bm_fit <- function(object, newdata = NULL, responses = NULL,
                           type = c("expected", "location", "predictive"),
                           random = "conditional",
                           new_levels = c("simulate", "zero", "error"),
                           spatial = c("conditional", "none"), n_draws = 1000,
                           level = 0.95, summarise = TRUE, seed = NULL, ...) {
  type <- match.arg(type)
  new_levels <- match.arg(new_levels)
  spatial <- match.arg(spatial)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  core <- bm_predict_core(object, newdata, responses, random, new_levels, spatial,
                          n_draws, seed)
  fams <- bm_families()
  draws <- stats::setNames(lapply(core$responses, function(r) {
    fam <- fams[[object$spec$family[[r]]]]
    switch(type,
           location = core$mu[[r]],
           expected = fam$expected(core$mu[[r]], core$sigma[[r]]),
           predictive = if (is.null(seed)) {
             fam$simulate(core$mu[[r]], core$sigma[[r]])
           } else {
             bm_with_seed(seed + 1L, fam$simulate(core$mu[[r]], core$sigma[[r]]))
           })
  }), core$responses)
  inp <- core$inputs
  n_na <- sum(!stats::complete.cases(t(draws[[1]])))
  if (n_na > 0L) {
    message(n_na, " row(s) have missing values in required inputs; predictions are NA.")
  }
  if (!summarise) {
    return(structure(list(draws = draws, data = inp$base, .row = inp$rows, type = type,
                          draw_index = core$draw_index, level = level),
                     class = "bm_draws"))
  }
  parts <- lapply(core$responses, function(r) {
    b <- inp$base
    b$.row <- inp$rows
    b$response <- r
    st <- bm_summarise_draws(draws[[r]], level)
    for (nm in names(st)) b[[nm]] <- st[[nm]]
    b
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  attr(out, "type") <- type
  attr(out, "observations") <- is.null(newdata)
  out
}
