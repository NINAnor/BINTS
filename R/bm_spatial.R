# Structured spatial effects.
#
# Spatial units are localities (one coordinate pair per locality, constant over
# years and replicates). Both approaches therefore define a field at the sites;
# predictions at new locations use the posterior of the field parameters.
#
# autocorrelation : nearest-neighbour Gaussian process (NNGP; Datta et al. 2016)
#   with exponential correlation exp(-d / ell). The sites are ordered by
#   longitude and each site conditions on its m nearest *earlier* sites, so the
#   joint density is a product of sparse conditionals
#       w_i | w_N(i) ~ N(B_i w_N(i), sd^2 F_i).
#   B_i and F_i depend on the length scale ell. To avoid linear algebra inside
#   the MCMC, ell is restricted to a log-spaced grid and B_i, F_i are
#   precomputed for each grid value; the index of the grid value is sampled
#   with a categorical sampler. Cost per sweep is O(L m), against O(L^3) for a
#   full GP.
#
# coordinates : thin-plate regression spline of the (planar, km) coordinates
#   (mgcv, bs = "tp"). The penalty null space (the linear plane) gets vague
#   fixed-effect priors; the penalised part is a random effect with one
#   smoothing sd, in the diagonal parameterisation used for the cyclic spline.

bm_is_lonlat <- function(crs) isTRUE(sf::st_is_longlat(sf::st_crs(crs)))

# Distances in km between rows of A and rows of B
bm_distance_km <- function(A, B, lonlat, km_per_unit = 1) {
  A <- matrix(A, ncol = 2)
  B <- matrix(B, ncol = 2)
  if (lonlat) {
    rad <- pi / 180
    la1 <- A[, 2] * rad
    la2 <- B[, 2] * rad
    a <- sin(outer(la1, la2, "-") / 2)^2 +
      outer(cos(la1), cos(la2)) * sin(outer(A[, 1] * rad, B[, 1] * rad, "-") / 2)^2
    6371 * 2 * asin(pmin(sqrt(a), 1))
  } else {
    sqrt(outer(A[, 1], B[, 1], "-")^2 + outer(A[, 2], B[, 2], "-")^2) * km_per_unit
  }
}

bm_resolve_spatial_options <- function(spatial_options) {
  opts <- list(m = 10L, n_range = 20L, range = NULL, k = 25L, crs = 4326,
               km_per_unit = 0.001)
  if (!is.null(spatial_options)) {
    if (!is.list(spatial_options) || !all(names(spatial_options) %in% names(opts))) {
      stop("`spatial_options` must be a list with elements from: ",
           paste(names(opts), collapse = ", "), ".", call. = FALSE)
    }
    opts[names(spatial_options)] <- spatial_options
  }
  int_ok <- function(x, min) is.numeric(x) && length(x) == 1L && x >= min &&
    x == round(x)
  if (!int_ok(opts$m, 1)) stop("`spatial_options$m` must be an integer >= 1.", call. = FALSE)
  if (!int_ok(opts$n_range, 2)) stop("`spatial_options$n_range` must be an integer >= 2.", call. = FALSE)
  if (!int_ok(opts$k, 5)) stop("`spatial_options$k` must be an integer >= 5.", call. = FALSE)
  if (!is.null(opts$range) && (!is.numeric(opts$range) || length(opts$range) != 2L ||
                               any(opts$range <= 0) || opts$range[1] >= opts$range[2])) {
    stop("`spatial_options$range` must be c(min, max) length scales in km.",
         call. = FALSE)
  }
  opts$lonlat <- bm_is_lonlat(opts$crs)
  opts
}

# Planar coordinates (km, centred, common scale) used by the spline
bm_plane_info <- function(xy, lonlat, km_per_unit) {
  info <- list(lonlat = lonlat, km_per_unit = km_per_unit,
               x0 = mean(xy[, 1]), y0 = mean(xy[, 2]), scale = 1)
  P <- bm_plane(xy, info)
  sc <- sqrt(mean(c(stats::var(P[, 1]), stats::var(P[, 2]))))
  if (!is.finite(sc) || sc == 0) {
    stop("Site coordinates have no spatial variation.", call. = FALSE)
  }
  info$scale <- sc
  info
}

bm_plane <- function(xy, info) {
  xy <- matrix(xy, ncol = 2)
  if (info$lonlat) {
    x <- (xy[, 1] - info$x0) * cos(info$y0 * pi / 180) * 111.32
    y <- (xy[, 2] - info$y0) * 110.57
  } else {
    x <- (xy[, 1] - info$x0) * info$km_per_unit
    y <- (xy[, 2] - info$y0) * info$km_per_unit
  }
  cbind(x, y) / info$scale
}

bm_tp_basis <- function(xy, opts) {
  info <- bm_plane_info(xy, opts$lonlat, opts$km_per_unit)
  P <- bm_plane(xy, info)
  k <- min(opts$k, nrow(unique(round(P, 8))) - 1L)
  if (k < 5L) {
    stop("At least 6 distinct site locations are needed for a coordinate ",
         "spline.", call. = FALSE)
  }
  sm <- mgcv::smoothCon(mgcv::s(cx, cy, bs = "tp", k = k),
                        data = data.frame(cx = P[, 1], cy = P[, 2]),
                        absorb.cons = TRUE)[[1]]
  eig <- eigen(sm$S[[1]], symmetric = TRUE)
  pos <- eig$values > max(eig$values) * 1e-8
  Tm <- cbind(
    eig$vectors[, pos, drop = FALSE] %*% diag(1 / sqrt(eig$values[pos]), nrow = sum(pos)),
    eig$vectors[, !pos, drop = FALSE]
  )
  list(basis = sm$X %*% Tm, K = ncol(Tm), Kp = sum(pos), smooth = sm,
       transform = Tm, info = info)
}

# Length-scale grid (km): from the median nearest-neighbour distance to half the
# extent of the sites, log-spaced.
bm_range_grid <- function(xy, opts) {
  if (!is.null(opts$range)) {
    return(exp(seq(log(opts$range[1]), log(opts$range[2]), length.out = opts$n_range)))
  }
  L <- nrow(xy)
  nn <- unlist(lapply(split(seq_len(L), ceiling(seq_len(L) / 500)), function(ii) {
    D <- bm_distance_km(xy[ii, , drop = FALSE], xy, opts$lonlat, opts$km_per_unit)
    D[cbind(seq_along(ii), ii)] <- Inf
    apply(D, 1, min)
  }))
  r <- apply(xy, 2, range)
  dmax <- max(
    bm_distance_km(r[c(1, 3)], r[c(2, 4)], opts$lonlat, opts$km_per_unit),
    bm_distance_km(r[c(1, 4)], r[c(2, 3)], opts$lonlat, opts$km_per_unit)
  )
  lo <- max(stats::median(nn[is.finite(nn)]), dmax * 0.005)
  hi <- max(dmax / 2, lo * 2)
  exp(seq(log(lo), log(hi), length.out = opts$n_range))
}

# Neighbour sets and per-grid-value conditional coefficients for the NNGP
bm_nngp_setup <- function(xy, opts) {
  L <- nrow(xy)
  if (L < 5L) stop("At least 5 sites are needed for spatial autocorrelation.", call. = FALSE)
  m <- min(opts$m, L - 1L)
  grid <- bm_range_grid(xy, opts)
  G <- length(grid)
  nbr <- matrix(1L, L, m)
  nb_info <- vector("list", L)
  for (i in 2:L) {
    d_i <- as.numeric(bm_distance_km(xy[i, , drop = FALSE], xy[seq_len(i - 1L), , drop = FALSE],
                                     opts$lonlat, opts$km_per_unit))
    nb <- order(d_i)[seq_len(min(m, i - 1L))]
    nbr[i, seq_along(nb)] <- nb
    nb_info[[i]] <- list(
      nb = nb, dvec = d_i[nb],
      Dnn = bm_distance_km(xy[nb, , drop = FALSE], xy[nb, , drop = FALSE],
                           opts$lonlat, opts$km_per_unit)
    )
  }
  coefs <- bm_nngp_coefs(nb_info, grid, m)
  # children[i]: the sites that have i among their neighbours
  par <- unlist(lapply(2:L, function(i) nb_info[[i]]$nb))
  kid <- rep(2:L, vapply(2:L, function(i) length(nb_info[[i]]$nb), integer(1)))
  nChild <- tabulate(par, nbins = L)
  maxc <- max(1L, max(nChild))
  child <- matrix(0L, L, maxc)
  for (p in seq_len(L)) if (nChild[p] > 0L) child[p, seq_len(nChild[p])] <- kid[par == p]
  storage.mode(nbr) <- "integer"
  list(m = m, G = G, range_grid = grid, nbr = nbr, Bnn = coefs$B,
       # B[l, j, g] flattened to a matrix with row (l - 1) m + j
       Bflat = matrix(aperm(coefs$B, c(2, 1, 3)), L * m, G),
       sqrtF = coefs$sqrtF, child = child, nChild = as.integer(nChild),
       xy = xy, options = opts)
}

bm_nngp_coefs <- function(nb_info, grid, m) {
  L <- length(nb_info)
  G <- length(grid)
  B <- array(0, c(L, m, G))
  sqrtF <- matrix(1, L, G)
  for (g in seq_len(G)) {
    for (i in 2:L) {
      ni <- nb_info[[i]]
      k <- length(ni$nb)
      C <- exp(-ni$Dnn / grid[g]) + diag(1e-6, k)
      cv <- exp(-ni$dvec / grid[g])
      b <- solve(C, cv)
      B[i, seq_len(k), g] <- b
      sqrtF[i, g] <- sqrt(max(1 - sum(cv * b), 1e-8))
    }
  }
  list(B = B, sqrtF = sqrtF)
}

bm_prepare_spatial <- function(dat, columns, spatial_effect, spatial_options) {
  if (spatial_effect == "none") return(NULL)
  opts <- bm_resolve_spatial_options(spatial_options)
  xcol <- columns$longitude
  ycol <- columns$latitude
  if (is.null(xcol) || is.null(ycol)) {
    stop("`columns` needs 'longitude' and 'latitude' entries (or supply ",
         "`coordinates`).", call. = FALSE)
  }
  if (!is.numeric(dat[[xcol]]) || !is.numeric(dat[[ycol]])) {
    stop("Coordinate columns '", xcol, "' and '", ycol, "' must be numeric.",
         call. = FALSE)
  }
  loc <- as.character(dat[[columns$locality]])
  xy_all <- cbind(dat[[xcol]], dat[[ycol]])
  first <- !duplicated(loc)
  site_names <- loc[first]
  xy <- xy_all[first, , drop = FALSE]
  mm <- match(loc, site_names)
  if (any(abs(xy_all - xy[mm, , drop = FALSE]) > 1e-8)) {
    stop("Coordinates must be constant within each locality.", call. = FALSE)
  }
  ord <- order(xy[, 1], xy[, 2])
  site_names <- site_names[ord]
  xy <- xy[ord, , drop = FALSE]
  out <- list(type = spatial_effect, site_names = site_names, xy = xy,
              site_id = match(loc, site_names), options = opts)
  if (spatial_effect == "autocorrelation") {
    out$nngp <- bm_nngp_setup(xy, opts)
  } else {
    out$tp <- bm_tp_basis(xy, opts)
  }
  out
}

# Model code for one spatial field of one response. `tag` names the field:
# "sp" is the main spatial effect; "spt" is the space-time interaction field
# (a spatially varying slope or seasonal amplitude), which reuses the same
# NNGP/spline machinery with its own parameters. `mult` is an optional
# expression the field is multiplied by in the linear predictor.
# Returns statements, the linear-predictor term, constants and monitors.
bm_spatial_block <- function(r, prep, scale, null_sd, tag = "sp", mult = NULL,
                              derived = NULL, share_to = NULL) {
  S <- function(p) as.name(paste(p, tag, r, sep = "_"))
  nm <- function(p) paste(p, tag, r, sep = "_")
  sp <- prep$spatial
  const <- stats::setNames(list(1 / scale^2), nm("prior_tau"))
  sd_stmt <- bquote(.(S("sd")) ~ T(dt(0, tau = .(S("prior_tau")), df = 3), 0, Inf))
  wrap <- function(field) {
    term <- bquote(.(field)[site_id[i]])
    if (is.null(mult)) term else bquote(.(term) * .(mult))
  }

  if (!is.null(derived)) {
    # Second response of a shared latent field: effect = lambda * reference field
    Sref <- function(p) as.name(paste(p, tag, derived$ref, sep = "_"))
    ld <- if (!is.null(derived$loading_sd)) bm_loading(tag, r, derived$loading_sd)
    lam <- function(e) if (is.null(ld)) e else bquote(.(ld$sym) * .(e))
    sd_expr <- if (is.null(ld)) Sref("sd") else bquote(abs(.(ld$sym)) * .(Sref("sd")))
    stmts <- c(if (!is.null(ld)) list(ld$stmt), list(bquote(.(S("sd")) <- .(sd_expr))))
    if (sp$type == "autocorrelation") {
      w <- S("w")
      # scalar nodes, so that each site's derived value has its own dependents
      stmts <- c(stmts, list(
        bquote(.(S("phi_idx")) <- .(Sref("phi_idx"))),
        bquote(for (l in 1:L_site) .(w)[l] <- .(lam(bquote(.(Sref("w"))[l]))))
      ))
      return(list(stmts = stmts, terms = list(wrap(w)),
                  constants = if (!is.null(ld)) ld$const else list(),
                  monitors = c(nm(c("sd", "phi_idx", "w")), ld$name),
                  field = nm("w")))
    }
    stmts <- c(stmts, list(
      bquote(.(S("v"))[1:K_sp] <- .(lam(bquote(.(Sref("v"))[1:K_sp])))),
      bquote(.(S("f"))[1:L_site] <- Bsp[1:L_site, 1:K_sp] %*% .(S("v"))[1:K_sp])
    ))
    return(list(stmts = stmts, terms = list(wrap(S("f"))),
                constants = if (!is.null(ld)) ld$const else list(),
                monitors = c(nm(c("sd", "v")), ld$name), field = nm("f")))
  }

  if (sp$type == "autocorrelation") {
    m <- sp$nngp$m
    w <- S("w")
    idx <- S("phi_idx")
    mean_expr <- Reduce(function(a, b) call("+", a, b), lapply(seq_len(m), function(j) {
      bquote(Bnn[l, .(j), .(idx)] * .(w)[nbr_idx[l, .(j)]])
    }))
    stmts <- list(
      sd_stmt,
      bquote(.(idx) ~ dcat(phi_prob[1:G_rng])),
      # one vector node with a custom density and its own sampler (see
      # bm_distributions.R)
      bquote(.(w)[1:L_site] ~ dnngp(sdev = .(S("sd")), idx = .(idx),
                                    Bflat = Bflat[1:LM_nn, 1:G_rng],
                                    sqrtF = sqrtF[1:L_site, 1:G_rng],
                                    nbr = nbr_idx[1:L_site, 1:m_nn]))
    )
    return(list(stmts = stmts, terms = list(wrap(w)), constants = const,
                monitors = nm(c("sd", "phi_idx", "w")), field = nm("w"),
                nngp = list(list(
                  target = nm("w"), sd = nm("sd"), idx = nm("phi_idx"),
                  mu = paste0("mu_", r),
                  derived = if (!is.null(share_to)) paste("w", tag, share_to, sep = "_"),
                  scale0 = 0.3 * scale))))
  }

  v <- S("v")
  f <- S("f")
  stmts <- list(
    sd_stmt,
    bquote(for (k in 1:Kp_sp) .(v)[k] ~ dnorm(0, sd = .(S("sd")))),
    bquote(for (k in (Kp_sp + 1):K_sp) .(v)[k] ~ dnorm(0, sd = .(S("prior_null_sd")))),
    bquote(.(f)[1:L_site] <- Bsp[1:L_site, 1:K_sp] %*% .(v)[1:K_sp])
  )
  const[[nm("prior_null_sd")]] <- null_sd
  list(stmts = stmts, terms = list(wrap(f)), constants = const,
       monitors = nm(c("sd", "v")), field = nm("f"))
}

# Constants shared by all responses
bm_spatial_constants <- function(prep) {
  sp <- prep$spatial
  out <- list(site_id = sp$site_id, L_site = nrow(sp$xy))
  if (sp$type == "autocorrelation") {
    ng <- sp$nngp
    out <- c(out, list(m_nn = ng$m, LM_nn = nrow(ng$xy) * ng$m, G_rng = ng$G, nbr_idx = ng$nbr,
                       phi_prob = rep(1 / ng$G, ng$G)))
  } else {
    out <- c(out, list(Bsp = sp$tp$basis, K_sp = sp$tp$K, Kp_sp = sp$tp$Kp))
  }
  out
}

# Posterior draws (iterations x nrow(new_xy)) of the spatial effect of
# `response` at new locations (columns: x = longitude/easting, y = latitude/
# northing, as in the fitting data). For the NNGP, `simulate = TRUE` draws the
# field at the new location from its conditional distribution given the
# fitted field (includes spatial uncertainty); FALSE returns the conditional
# mean. The coordinate spline is deterministic given the coefficients.
bm_spatial_draws <- function(fit, response, new_xy, simulate = TRUE, tag = "sp",
                             mat = NULL) {
  sp <- fit$prep$spatial
  if (is.null(sp)) stop("The model has no spatial effect.", call. = FALSE)
  new_xy <- matrix(as.matrix(new_xy), ncol = 2)
  if (is.null(mat)) mat <- as.matrix(fit$samples)
  if (sp$type == "coordinates") {
    tp <- sp$tp
    P <- bm_plane(new_xy, tp$info)
    B <- mgcv::PredictMat(tp$smooth, data.frame(cx = P[, 1], cy = P[, 2])) %*% tp$transform
    v <- mat[, paste0("v_", tag, "_", response, "[", seq_len(tp$K), "]"), drop = FALSE]
    return(v %*% t(B))
  }
  ng <- sp$nngp
  opts <- ng$options
  L <- nrow(ng$xy)
  W <- mat[, paste0("w_", tag, "_", response, "[", seq_len(L), "]"), drop = FALSE]
  gidx <- mat[, paste0("phi_idx_", tag, "_", response)]
  sdv <- mat[, paste0("sd_", tag, "_", response)]
  out <- matrix(0, nrow(mat), nrow(new_xy))
  for (p in seq_len(nrow(new_xy))) {
    d <- as.numeric(bm_distance_km(new_xy[p, , drop = FALSE], ng$xy,
                                   opts$lonlat, opts$km_per_unit))
    nb <- order(d)[seq_len(ng$m)]
    Dnn <- bm_distance_km(ng$xy[nb, , drop = FALSE], ng$xy[nb, , drop = FALSE],
                          opts$lonlat, opts$km_per_unit)
    for (g in unique(gidx)) {
      rows <- which(gidx == g)
      ell <- ng$range_grid[g]
      cv <- exp(-d[nb] / ell)
      b <- solve(exp(-Dnn / ell) + diag(1e-6, length(nb)), cv)
      f <- max(1 - sum(cv * b), 1e-8)
      mu <- W[rows, nb, drop = FALSE] %*% b
      out[rows, p] <- mu + if (simulate) stats::rnorm(length(rows), 0, sdv[rows] * sqrt(f)) else 0
    }
  }
  out
}
