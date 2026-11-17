# Prior specification. Priors are weakly informative and set relative to the
# scale of each response (the response itself for the truncated normal, its
# logarithm for the log-normal). Covariates are standardised internally, so
# `beta_sd` is the prior sd of an effect per 1 sd of the covariate.
#
#   beta0      ~ Normal(beta0_mean, beta0_sd)
#   beta[k]    ~ Normal(0, beta_sd)
#   sd_re      ~ half-t(3, 0, sd_re_scale)    (random-effect sds)
#   sigma      ~ half-t(3, 0, sigma_scale)    (residual sd)
#   trend      ~ Normal(0, trend_sd)          (linear trend, per year)
#   sd_cyc     ~ half-t(3, 0, smooth_scale)   (cyclic spline smoothing sd)
#   sd_sp      ~ half-t(3, 0, spatial_scale)  (NNGP marginal sd / spatial spline sd)
#   sd_spt     ~ half-t(3, 0, interaction_scale)  (space-time interaction field: sd of
#                 spatially varying slope per year, or of the amplitude multiplier)
#   phi        ~ Gamma(disp_shape, disp_rate) (negative binomial size; mean 10 by default)
#   lambda     ~ Normal(0, loading_sd)        (loading of the second response on a shared
#                 latent effect; loading_sd = 5 x response-scale ratio)
#   range      ~ uniform over a log-spaced grid of length scales (NNGP)

bm_prior_names <- function() {
  c("beta0_mean", "beta0_sd", "beta_sd", "sd_re_scale", "sigma_scale",
    "trend_sd", "smooth_scale", "spatial_scale", "interaction_scale",
    "loading_sd", "disp_shape", "disp_rate")
}

bm_default_priors <- function(prep, family) {
  fam <- bm_families()
  scales <- vapply(names(prep$y), function(r) {
    stats::sd(fam[[family[[r]]]]$prior_scale(prep$y[[r]]))
  }, numeric(1))
  lapply(stats::setNames(names(prep$y), names(prep$y)), function(r) {
    z <- fam[[family[[r]]]]$prior_scale(prep$y[[r]])
    s <- stats::sd(z)
    if (!is.finite(s) || s == 0) {
      stop("Response '", r, "' has no variation.", call. = FALSE)
    }
    yr_sd <- if (identical(prep$temporal$type, "linear")) prep$temporal$year_sd else 1
    list(beta0_mean = mean(z), beta0_sd = 10 * s, beta_sd = 2.5 * s,
         sd_re_scale = s, sigma_scale = s,
         trend_sd = 2.5 * s / yr_sd, smooth_scale = s,
         spatial_scale = s,
         loading_sd = 5 * s / scales[[1]],
         disp_shape = 1, disp_rate = 0.1,
         interaction_scale = if (identical(prep$temporal$type, "linear")) 2.5 * s / yr_sd else 1)
  })
}

# Merge user priors (list keyed by response, each a named list of overrides)
# with defaults.
bm_resolve_priors <- function(priors, defaults) {
  if (is.null(priors)) return(defaults)
  if (!is.list(priors) || is.null(names(priors)) ||
      !all(names(priors) %in% names(defaults))) {
    stop("`priors` must be a list named by response (",
         paste(names(defaults), collapse = ", "), "), each a named list with ",
         "elements from: ", paste(bm_prior_names(), collapse = ", "), ".",
         call. = FALSE)
  }
  for (r in names(priors)) {
    bad <- setdiff(names(priors[[r]]), bm_prior_names())
    if (length(bad)) {
      stop("Unknown prior element(s) for '", r, "': ", paste(bad, collapse = ", "),
           ". Allowed: ", paste(bm_prior_names(), collapse = ", "), ".",
           call. = FALSE)
    }
    vals <- unlist(priors[[r]])
    if (!is.numeric(vals) || anyNA(vals) ||
        any(vals[names(vals) != "beta0_mean"] <= 0)) {
      stop("Prior values for '", r, "' must be numeric, and scales positive.",
           call. = FALSE)
    }
    defaults[[r]][names(priors[[r]])] <- priors[[r]]
  }
  defaults
}
