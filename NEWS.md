# BINTS 0.2.0

* New `bm_*` modelling framework built on NIMBLE: `fit_model()` fits joint
  biomass and richness models.
* Likelihoods: truncated normal or lognormal for biomass; truncated normal,
  Poisson or negative binomial for richness.
* Spatial effects: nearest-neighbour Gaussian process (NNGP) or thin-plate
  spline. Temporal effects: linear trend or penalised cyclic cubic spline,
  with an optional space-time interaction.
* Parameters, latent effects and random effects can be shared between the
  responses, with free loadings or identical values.
* `predict()`, `posterior_predictive_check()`, `residuals()` and plotting
  helpers for fitted models; `simulate_bm_data()` for simulation.
* `scripts/allInsects_fit.R` fits the all-insects joint model (truncated normal
  biomass, negative binomial richness, NNGP spatial field, cyclic seasonal
  spline, shared effects) and produces posterior predictive checks, residual
  plots and prediction maps. `scripts/run_allInsects_server.R` runs it from the
  command line.
* The workflow vignette now uses a small simulated data set so that it renders
  in a few minutes.
