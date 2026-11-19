# BINTS

Bayesian hierarchical spatio-temporal models of insect biomass and species
richness, fitted with [NIMBLE](https://r-nimble.org). Biomass and richness can
be modelled separately or jointly, with effects shared between the two.

## Installation

NIMBLE compiles C++ at run time, so a working C++ toolchain is required
(Rtools on Windows).

```r
# install.packages("pak")
pak::pak("NINAnor/BINTS")      # or devtools::install() from a local clone

# For the development version:
pak::pak("NINAnor/BINTS@devel")  
```

## Model

For response $r$ (biomass or richness), observation $i$ in year $t$, locality
$l$ and replicate $k$:

$$\mu_{ri} = \beta_{0r} + x_i^\top\beta_r + f_r(\text{time}) + s_r(\text{space}) +
u^{year}_{rt} + u^{locality}_{rl} + u^{rep}_{rk}$$

| Component | Options |
|---|---|
| Biomass likelihood | truncated normal (truncation at zero in the likelihood), lognormal |
| Richness likelihood | negative binomial, Poisson, truncated normal |
| Spatial effect | nearest-neighbour Gaussian process (NNGP), thin-plate spline of coordinates |
| Temporal effect | linear trend, penalised cyclic cubic spline (seasonality) |
| Space-time | optional spatially varying slope or seasonal amplitude |
| Random effects | year, locality, replication |
| Joint models | shared parameters, shared latent processes with loadings, correlated random effects |

Sampling uses NIMBLE MCMC with recorded seeds, so fits are reproducible.

## Quick start

```r
library(BINTS)

sim <- simulate_bm_data(n_year = 4, n_locality = 20, n_rep = 2, seed = 1)

fit <- fit_model(
  sim,
  fixed_effects = c("x1", "x2"),
  temporal_effect = "linear",
  spatial_effect = "coordinates",
  n_chains = 2, n_iter = 2000, n_burnin = 1000, n_thin = 2,
  parallel = FALSE
)

summary(fit)
posterior_predictive_check(fit)
predict(fit)
```

For real data with counts of species, use
`family = c(biomass = "truncnorm", richness = "negbinomial")`, and for an NNGP
field use `spatial_effect = "autocorrelation"`. Data need one row per sampling
event with the columns listed by `bm_default_columns()` (override with
`columns =`).

## Main functions

| Function | Purpose |
|---|---|
| `fit_model()` | Fit a biomass and/or richness model |
| `predict()` | Posterior predictions at observed or new locations |
| `posterior_predictive_check()` | Posterior predictive checks and discrepancy summaries |
| `residuals()` | Raw, standardised and posterior-predictive residuals |
| `plot_observed_predicted()`, `plot_residuals()`, `plot_response_covariates()`, `plot_spatial_predictions()` | Diagnostic plots and maps |
| `simulate_bm_data()` | Simulate data with known parameters |

The earlier `create_insect_*()`, `run_insect_nimble_model()` and
`select_insect_module()` functions are retained.

## Practical notes

* Compilation and sampling are slow. Defaults (4 chains, 50000 iterations) are
  for real analyses; use fewer iterations for testing.
* NNGP cost grows quickly with the number of localities.
* A cyclic seasonal effect and a replication random effect can be partly
  confounded when replicates correspond to sampling occasions.

## Documentation

See `vignette("BINTS-workflow")` for a full workflow, `?fit_model` for all
options, and [NEWS.md](NEWS.md) for changes.

## Licence

MIT
