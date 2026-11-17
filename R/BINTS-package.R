#' BINTS: Bayesian hierarchical spatio-temporal models of insect biomass and
#' species richness
#'
#' Tools to fit, check, predict from and visualise Bayesian hierarchical models
#' of insect biomass and species richness with NIMBLE.
#'
#' @section Getting started:
#' ```
#' fit <- fit_model(my_data,
#'   responses = c("biomass", "richness"),
#'   fixed_effects = c("mean_temp", "mean_precip"),
#'   temporal_effect = "linear", spatial_effect = "coordinates",
#'   seed = 1)
#' fit
#' summary(fit)
#' posterior_predictive_check(fit)
#' predict(fit, newdata)
#' ```
#' Column names are mapped through [bm_default_columns()] (for example
#' `sum_wet_weight` for biomass and `no_species` for richness); pass a modified
#' list to `fit_model(columns = )` for other names. The vignette
#' `vignette("BINTS-workflow", package = "BINTS")` walks through a complete
#' analysis on simulated data ([simulate_bm_data()]).
#'
#' @section The model:
#' For response \eqn{r} (biomass or richness) and observation \eqn{i} in year
#' \eqn{t}, locality \eqn{l}, replicate \eqn{k},
#' \deqn{y_{ri} \sim \mathrm{TruncNormal}(\mu_{ri}, \sigma_r^2;\ 0, \infty),}
#' \deqn{\mu_{ri} = \beta_{0r} + x_i^\top \beta_r + f_r(\text{time}) +
#'   s_r(\text{space}) + g_r(\text{space}, \text{time}) + u^{year}_{rt} +
#'   u^{locality}_{rl} + u^{rep}_{rk}.}
#' The truncation at zero is part of the likelihood, so predictions are never
#' negative. `mu` and `sigma` are the location and scale of the *untruncated*
#' normal, not the mean and sd of the response; [predict.bm_fit()] returns the
#' mean of the truncated distribution (the *expected response*). Richness is a
#' count: the truncated normal ignores its discreteness and is adequate for
#' large counts. The family registry is designed so that count families can be
#' added.
#'
#' @section Main functions:
#' | Task | Functions |
#' |:-----|:----------|
#' | Fit | [fit_model()], with [bm_default_columns()] for column names |
#' | Inspect | `print()`, `summary()` ([bm_fit-methods]) |
#' | Predict | [predict.bm_fit()] (expected response, location, predictive draws; conditional or marginal over random effects; new years, locations, grids) |
#' | Check | [posterior_predictive_check()], [residuals.bm_fit()] |
#' | Map | [plot_spatial_predictions()] |
#' | Diagnose | [plot_observed_predicted()], [plot_residuals()], [plot_response_covariates()] |
#' | Simulate | [simulate_bm_data()] |
#'
#' @section Temporal, spatial and joint structures:
#' See the sections of [fit_model()]: *Temporal effects* (linear trend or cyclic
#' spline, and how a structured trend differs from the unstructured year
#' effect), *Spatial effects* (nearest-neighbour Gaussian process or coordinate
#' spline, their assumptions and requirements), *Spatio-temporal models*
#' (additive, or with a space-time interaction) and *Joint models* (sharing a
#' parameter, sharing a latent process with a loading, or correlated effects).
#'
#' @section Reproducibility:
#' Every chain gets its own seed, derived from `seed` and saved in `fit$seeds`;
#' supply `chain_seeds` to reuse them. The same data, specification and seeds
#' give identical samples whether chains run in parallel or not. Fitted objects
#' also store the priors, MCMC settings, model code and monitored parameters.
#'
#' @section Earlier functions:
#' The `create_insect_*()` functions, [select_insect_module()],
#' [run_insect_nimble_model()] and [prepare_cubic_cyclic_temporal_spline()]
#' from earlier versions are unchanged and still exported.
#'
#' @keywords internal
#' @importFrom stats predict residuals na.omit
#' @importFrom ggplot2 .data
#' @md
"_PACKAGE"
