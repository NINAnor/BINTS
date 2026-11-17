# Posterior predictive checks ---------------------------------------------------

bm_rowwise_stats <- function(Y, quantiles, near_zero) {
  m <- ncol(Y)
  mu <- rowMeans(Y)
  out <- cbind(
    mean = mu,
    sd = if (m > 1L) sqrt(pmax(rowSums((Y - mu)^2) / (m - 1L), 0)) else NA_real_,
    min = apply(Y, 1, min),
    max = apply(Y, 1, max)
  )
  q <- t(apply(Y, 1, stats::quantile, probs = quantiles, names = FALSE))
  if (length(quantiles) == 1L) q <- matrix(q, ncol = 1L)
  colnames(q) <- paste0("q", format(100 * quantiles, trim = TRUE, drop0trailing = TRUE))
  cbind(out, q, prop_near_zero = rowMeans(Y <= near_zero), prop_zero = rowMeans(Y == 0))
}

bm_midp <- function(t_rep, t_obs) {
  if (all(is.na(t_rep)) || all(is.na(t_obs))) return(NA_real_)
  mean(t_rep > t_obs) + 0.5 * mean(t_rep == t_obs)
}

#' Posterior predictive checks for a fitted model
#'
#' Simulates replicated data from the fitted posterior predictive distribution
#' and compares discrepancy statistics of the replicates with those of the
#' observed data, overall and optionally within groups.
#'
#' @section Procedure:
#' For each retained posterior draw \eqn{\theta} (all parameters, including the
#' estimated random effects and spatial fields) the linear predictor
#' \eqn{\mu_{ri}(\theta)} is computed at the observed rows and a replicate
#' \eqn{y^{rep}_{ri}} is drawn from **the same likelihood and truncation as the
#' fitted model** (truncated normal with lower bound 0, or the log-normal), so
#' replicates are never negative. Posterior uncertainty is propagated because
#' every replicate uses its own draw \eqn{\theta}. With `random =
#' "marginal"` the random effects are redrawn from their population
#' distribution instead of being conditioned on, which checks the model for new
#' groups (a mixed predictive check).
#'
#' This is a *posterior predictive observation*, not the posterior of the
#' expected response: `predict(type = "expected")` returns the latter, whose
#' uncertainty excludes observation noise. Predictive replicates include it.
#'
#' @section Discrepancies:
#' Data statistics (computed on observed and on replicated data): `mean`, `sd`,
#' `min`, `max`, the requested `quantiles` (`q5`, `q50`, ...), the proportion of
#' values at or below `near_zero` (`prop_near_zero`) and the proportion of exact
#' zeros (`prop_zero`; always 0 under the truncated normal, so observed zeros
#' are flagged). Likelihood-based discrepancies depend on the parameters and
#' are evaluated as \eqn{T(y, \theta)} for the data and \eqn{T(y^{rep}, \theta)}
#' for the replicates with the same draw: `"deviance"`,
#' \eqn{-2\sum_i \log p(y_i \mid \theta)}, and `"chisq"`, the Pearson-type
#' \eqn{\sum_i (y_i - E_i)^2 / V_i} with the mean and variance of the truncated
#' distribution.
#'
#' The Bayesian posterior predictive p-value is
#' \eqn{\Pr(T^{rep} \ge T^{obs} \mid y)}, estimated over the draws with ties
#' counted half (mid-p), so it is well behaved for discrete statistics. Values
#' near 0 or 1 (flagged `"extreme"` outside 0.025-0.975) indicate that the model
#' does not reproduce that feature of the data. These p-values are not
#' uniformly distributed even for a correct model (they are conservative for
#' statistics the model is fitted to, such as the mean), so treat them as
#' diagnostics, not tests.
#'
#' @section Richness:
#' Richness is a count; the truncated normal generates continuous values with
#' no mass at zero. `prop_near_zero`, `prop_zero` and the lower quantiles
#' expose this when counts are small. Replicates are not rounded.
#'
#' @param model A `bm_fit` object.
#' @param responses Responses to check (default all fitted).
#' @param n_draws Number of posterior draws (evenly spaced) used for the check.
#' @param discrepancy `"deviance"` and/or `"chisq"`.
#' @param by `NULL`, or the names of one or more columns of the modelling data
#'   (for example `"year"`, `"locality"`, `"region_name"`); statistics are
#'   additionally computed within each group (for several columns, the
#'   combinations).
#' @param quantiles Quantiles of the data to compare.
#' @param near_zero Threshold for `prop_near_zero`; default 5 percent of the
#'   observed mean of each response.
#' @param random Random-effect mode passed to the prediction engine
#'   (`"conditional"` by default; see [predict.bm_fit()]).
#' @param n_overlay Number of replicated densities drawn in the overlay plot.
#' @param seed Seed for the replicates (the RNG state of the session is restored).
#'
#' @return An object of class `bm_ppc`, a list with `summary` (one row per
#'   response, group and statistic: `observed`, replicated mean and 95%
#'   interval, `p_value`, `flag`), `plots` (ggplot2 objects: `density`,
#'   `statistics`, `discrepancy`, and `groups` if `by` is given), `draws` (the
#'   discrepancy draws), `replicates` (the first `n_overlay` replicated data
#'   sets per response) and `settings`. `plot()` returns or draws the plots.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 12, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", fixed_effects = "x1",
#'                  n_chains = 2, n_iter = 2000, n_burnin = 1000, n_thin = 2,
#'                  seed = 1, parallel = FALSE)
#' ppc <- posterior_predictive_check(fit, n_draws = 200, by = "year", seed = 1)
#' ppc
#' ppc$plots$density
#' }
#' @md
posterior_predictive_check <- function(model,
                                       responses = NULL,
                                       n_draws = 1000,
                                       discrepancy = "deviance",
                                       by = NULL,
                                       quantiles = c(0.05, 0.5, 0.95),
                                       near_zero = NULL,
                                       random = "conditional",
                                       n_overlay = 50,
                                       seed = NULL) {
  if (!inherits(model, "bm_fit")) stop("`model` must be a bm_fit object.", call. = FALSE)
  if (!is.character(discrepancy) || !all(discrepancy %in% c("deviance", "chisq"))) {
    stop("`discrepancy` must be 'deviance' and/or 'chisq'.", call. = FALSE)
  }
  if (!is.numeric(quantiles) || any(quantiles <= 0 | quantiles >= 1)) {
    stop("`quantiles` must be numbers in (0, 1).", call. = FALSE)
  }
  prep <- model$prep
  dat <- prep$data
  if (!is.null(by)) {
    if (!is.character(by) || !all(by %in% names(dat))) {
      stop("`by` must name column(s) of the modelling data: ",
           paste(setdiff(by, names(dat)), collapse = ", "), " not found.", call. = FALSE)
    }
    grp <- as.character(interaction(dat[by], drop = TRUE, sep = " / "))
  } else grp <- NULL

  core <- bm_predict_core(model, NULL, responses, random, "simulate", "conditional",
                          n_draws, seed)
  resp <- core$responses
  fams <- bm_families()
  D <- nrow(core$mu[[1]])
  stat_names <- NULL
  rows <- list()
  dens <- list()
  scat <- list()
  hist_df <- list()
  reps_keep <- list()
  draws_out <- list()

  one <- function() {
    for (r in resp) {
      fam <- fams[[model$spec$family[[r]]]]
      y <- prep$y[[r]]
      n <- length(y)
      mu <- core$mu[[r]]
      sg <- core$sigma[[r]]
      Yrep <- fam$simulate(mu, sg)
      nz <- if (is.null(near_zero)) 0.05 * mean(y) else near_zero
      Lobs <- fam$loglik(matrix(y, D, n, byrow = TRUE), mu, sg)
      Lrep <- fam$loglik(Yrep, mu, sg)
      if ("chisq" %in% discrepancy) {
        Ev <- fam$expected(mu, sg)
        Vv <- fam$variance(mu, sg)
        Cobs <- (matrix(y, D, n, byrow = TRUE) - Ev)^2 / Vv
        Crep <- (Yrep - Ev)^2 / Vv
      }
      groups <- c(list(all = seq_len(n)),
                  if (!is.null(grp)) split(seq_len(n), grp))
      for (g in names(groups)) {
        ii <- groups[[g]]
        if (length(ii) < 2L && g != "all") next
        Srep <- bm_rowwise_stats(Yrep[, ii, drop = FALSE], quantiles, nz)
        Sobs <- bm_rowwise_stats(matrix(y[ii], 1L), quantiles, nz)
        T_obs <- matrix(Sobs, D, ncol(Sobs), byrow = TRUE, dimnames = list(NULL, colnames(Sobs)))
        T_rep <- Srep
        if ("deviance" %in% discrepancy) {
          T_obs <- cbind(T_obs, deviance = -2 * rowSums(Lobs[, ii, drop = FALSE]))
          T_rep <- cbind(T_rep, deviance = -2 * rowSums(Lrep[, ii, drop = FALSE]))
        }
        if ("chisq" %in% discrepancy) {
          T_obs <- cbind(T_obs, chisq = rowSums(Cobs[, ii, drop = FALSE]))
          T_rep <- cbind(T_rep, chisq = rowSums(Crep[, ii, drop = FALSE]))
        }
        for (s in colnames(T_rep)) {
          p <- bm_midp(T_rep[, s], T_obs[, s])
          rows[[length(rows) + 1L]] <<- data.frame(
            response = r, group = g, statistic = s,
            observed = mean(T_obs[, s]),
            replicated_mean = mean(T_rep[, s]),
            replicated_lower = stats::quantile(T_rep[, s], 0.025, names = FALSE, na.rm = TRUE),
            replicated_upper = stats::quantile(T_rep[, s], 0.975, names = FALSE, na.rm = TRUE),
            p_value = p,
            flag = if (!is.na(p) && (p < 0.025 || p > 0.975)) "extreme" else "",
            stringsAsFactors = FALSE
          )
        }
        if (g == "all") {
          draws_out[[r]] <<- list(observed = T_obs, replicated = T_rep)
          long <- data.frame(
            response = r,
            statistic = rep(colnames(T_rep), each = D),
            value = as.vector(T_rep),
            observed = as.vector(T_obs)
          )
          is_disc <- long$statistic %in% c("deviance", "chisq")
          hist_df[[r]] <<- long[!is_disc, , drop = FALSE]
          scat[[r]] <<- long[is_disc, , drop = FALSE]
        }
      }
      # Overlay densities of observed and replicated data
      k <- min(n_overlay, D)
      hi <- as.numeric(stats::quantile(c(y, Yrep[seq_len(k), ]), 0.995))
      xs <- seq(0, hi, length.out = 256)
      dn <- function(v) stats::density(v, from = 0, to = hi, n = 256)$y
      dens[[r]] <<- rbind(
        data.frame(response = r, x = xs, density = dn(y), series = "observed", id = 0L),
        do.call(rbind, lapply(seq_len(k), function(j) {
          data.frame(response = r, x = xs, density = dn(Yrep[j, ]), series = "replicated", id = j)
        }))
      )
      reps_keep[[r]] <<- Yrep[seq_len(k), , drop = FALSE]
    }
  }
  if (is.null(seed)) one() else bm_with_seed(seed + 2L, one())

  summ <- do.call(rbind, rows)
  rownames(summ) <- NULL
  plots <- bm_ppc_plots(summ, dens, hist_df, scat, !is.null(by))
  structure(
    list(summary = summ, plots = plots, draws = draws_out, replicates = reps_keep,
         settings = list(responses = resp, n_draws = D, discrepancy = discrepancy,
                         by = by, quantiles = quantiles, random = random,
                         draw_index = core$draw_index)),
    class = "bm_ppc"
  )
}

bm_ppc_plots <- function(summ, dens, hist_df, scat, has_groups) {
  dens <- do.call(rbind, dens)
  dens$line <- paste(dens$response, dens$series, dens$id)
  p_dens <- ggplot2::ggplot(dens, ggplot2::aes(x = .data$x, y = .data$density, group = .data$line,
                                               colour = .data$series, alpha = .data$series,
                                               linewidth = .data$series)) +
    ggplot2::geom_line() +
    ggplot2::scale_alpha_manual(values = c(observed = 1, replicated = 0.25)) +
    ggplot2::scale_linewidth_manual(values = c(observed = 1, replicated = 0.3)) +
    ggplot2::scale_colour_manual(values = c(observed = "black", replicated = "steelblue")) +
    ggplot2::facet_wrap(~ response, scales = "free") +
    ggplot2::labs(x = "Value", y = "Density", colour = NULL, alpha = NULL, linewidth = NULL,
                  title = "Observed data and replicated data sets")

  h <- do.call(rbind, hist_df)
  obs_lines <- unique(h[, c("response", "statistic", "observed")])
  p_stat <- ggplot2::ggplot(h, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    ggplot2::geom_vline(data = obs_lines, ggplot2::aes(xintercept = .data$observed),
                        colour = "firebrick", linewidth = 0.8) +
    ggplot2::facet_grid(response ~ statistic, scales = "free") +
    ggplot2::labs(x = "Replicated value (red line: observed)", y = "Draws",
                  title = "Data statistics")

  sc <- do.call(rbind, scat)
  p_disc <- NULL
  if (!is.null(sc) && nrow(sc)) {
    pv <- summ[summ$group == "all" & summ$statistic %in% unique(sc$statistic),
               c("response", "statistic", "p_value")]
    pv$label <- sprintf("p = %.2f", pv$p_value)
    p_disc <- ggplot2::ggplot(sc, ggplot2::aes(x = .data$observed, y = .data$value)) +
      ggplot2::geom_point(alpha = 0.3, size = 0.8) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "firebrick") +
      ggplot2::geom_text(data = pv, ggplot2::aes(label = .data$label), x = -Inf, y = Inf,
                         hjust = -0.1, vjust = 1.5, inherit.aes = FALSE) +
      ggplot2::facet_grid(response ~ statistic, scales = "free") +
      ggplot2::labs(x = "Discrepancy of observed data", y = "Discrepancy of replicated data",
                    title = "Likelihood-based discrepancies")
  }
  out <- list(density = p_dens, statistics = p_stat, discrepancy = p_disc)
  if (has_groups) {
    g <- summ[summ$group != "all" & summ$statistic %in% c("mean", "sd"), ]
    if (nrow(g)) {
      out$groups <- ggplot2::ggplot(g, ggplot2::aes(x = .data$group)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = .data$replicated_lower,
                                             ymax = .data$replicated_upper),
                                colour = "steelblue", linewidth = 1.5, alpha = 0.6) +
        ggplot2::geom_point(ggplot2::aes(y = .data$observed), colour = "black", size = 1.2) +
        ggplot2::facet_grid(statistic ~ response, scales = "free_y") +
        ggplot2::labs(x = NULL, y = "Value",
                      title = "Observed (points) versus replicated 95% intervals by group") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
  }
  Filter(Negate(is.null), out)
}

#' @rdname posterior_predictive_check
#' @param x A `bm_ppc` object.
#' @param which Which plot: one of `names(x$plots)` (default `"density"`).
#' @param ... Unused.
#' @export
#' @md
plot.bm_ppc <- function(x, which = "density", ...) {
  if (!which %in% names(x$plots)) {
    stop("`which` must be one of: ", paste(names(x$plots), collapse = ", "), ".",
         call. = FALSE)
  }
  x$plots[[which]]
}

#' @rdname posterior_predictive_check
#' @export
#' @md
print.bm_ppc <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("<bm_ppc> posterior predictive check\n")
  cat("  Responses :", paste(s$responses, collapse = ", "), "\n")
  cat("  Draws     :", s$n_draws, "| random effects:", s$random,
      if (!is.null(s$by)) paste0("| by: ", paste(s$by, collapse = ", ")), "\n")
  tab <- x$summary[x$summary$group == "all", ]
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = digits)
  print(tab[, c("response", "statistic", "observed", "replicated_mean",
                "replicated_lower", "replicated_upper", "p_value", "flag")],
        row.names = FALSE)
  if (!is.null(s$by)) {
    g <- x$summary[x$summary$group != "all", ]
    cat("\n  Group-level results:", length(unique(g$group)), "groups;",
        sum(g$flag == "extreme"), "of", nrow(g), "statistics flagged extreme.\n")
  }
  invisible(x)
}
