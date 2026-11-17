# Structured temporal effects.
#
# linear : trend_r * (year - mean(year)); slope is per year. Year is centred so
#          the trend is orthogonal to the intercept.
# cyclic : penalised cubic cyclic regression spline (mgcv, bs = "cc") of a
#          seasonal variable (default: day of year), with a sum-to-zero
#          constraint so it is not confounded with the intercept.
#
# The cyclic spline is written in a diagonal form (Wood 2017): with penalty
# S = U D U', f(t) = B U D^{-1/2} w, w ~ N(0, sd^2 I). The single smoothing sd
# has a half-t prior, which is the Bayesian analogue of the penalty. The
# centred form mixes better than a non-centred one when data are informative.

bm_cyclic_basis <- function(season, k, period) {
  if (any(season < 0 | season > period)) {
    stop("Cyclic variable must lie within [0, period = ", period, "]; found ",
         "range [", min(season), ", ", max(season), "]. Adjust ",
         "`temporal_options$period`.", call. = FALSE)
  }
  sm <- mgcv::smoothCon(
    mgcv::s(season, bs = "cc", k = k), data = data.frame(season = season),
    knots = list(season = c(0, period)), absorb.cons = TRUE
  )[[1]]
  eig <- eigen(sm$S[[1]], symmetric = TRUE)
  pos <- eig$values > max(eig$values) * 1e-8
  Tm <- eig$vectors[, pos, drop = FALSE] %*% diag(1 / sqrt(eig$values[pos]),
                                                 nrow = sum(pos))
  list(basis = sm$X %*% Tm, K = ncol(Tm), smooth = sm, transform = Tm,
       period = period)
}

bm_prepare_temporal <- function(dat, columns, temporal_effect, temporal_options) {
  if (temporal_effect == "none") return(NULL)
  if (temporal_effect == "linear") {
    yr <- dat[[columns$year]]
    if (!is.numeric(yr)) {
      stop("A linear temporal effect needs a numeric year column ('",
           columns$year, "').", call. = FALSE)
    }
    if (length(unique(yr)) < 2L) {
      stop("A linear temporal effect needs at least two distinct years.",
           call. = FALSE)
    }
    return(list(type = "linear", year_center = mean(yr),
                year_c = yr - mean(yr), year_sd = stats::sd(yr)))
  }
  sn <- dat[[columns$season]]
  if (!is.numeric(sn)) {
    stop("A cyclic temporal effect needs a numeric seasonal column ('",
         columns$season, "').", call. = FALSE)
  }
  cb <- bm_cyclic_basis(sn, temporal_options$k, temporal_options$period)
  list(type = "cyclic", cyclic = cb, season = sn)
}

bm_resolve_temporal_options <- function(temporal_options) {
  opts <- list(k = 8L, period = 365)
  if (!is.null(temporal_options)) {
    if (!is.list(temporal_options) ||
        !all(names(temporal_options) %in% names(opts))) {
      stop("`temporal_options` must be a list with elements from: ",
           paste(names(opts), collapse = ", "), ".", call. = FALSE)
    }
    opts[names(temporal_options)] <- temporal_options
  }
  if (!is.numeric(opts$k) || opts$k < 4 || opts$k != round(opts$k)) {
    stop("`temporal_options$k` must be an integer >= 4.", call. = FALSE)
  }
  if (!is.numeric(opts$period) || opts$period <= 0) {
    stop("`temporal_options$period` must be a positive number.", call. = FALSE)
  }
  opts
}

# Posterior draws (iterations x length(season)) of the cyclic effect of
# response `response` at the given seasonal values.
bm_cyclic_draws <- function(fit, response, season, mat = NULL) {
  cb <- fit$prep$temporal$cyclic
  if (is.null(cb)) stop("The model has no cyclic temporal effect.", call. = FALSE)
  B <- mgcv::PredictMat(cb$smooth, data.frame(season = season)) %*% cb$transform
  if (is.null(mat)) mat <- as.matrix(fit$samples)
  w <- mat[, paste0("w_cyc_", response, "[", seq_len(cb$K), "]"), drop = FALSE]
  w %*% t(B)
}
