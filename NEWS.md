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
* The workflow vignette now uses a small simulated data set so that it renders
  in a few minutes.
