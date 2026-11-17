# BINTS

Bayesian joint models for insect biomass and species richness with spatial
and temporal structure, fitted with NIMBLE.

## Installation

```r
# install.packages("pak")
pak::pak("<owner>/BINTS")   # or devtools::install() from a local clone
```

## Model

* Biomass: truncated normal or lognormal.
* Richness (a count): negative binomial, Poisson or truncated normal.
* Space: nearest-neighbour Gaussian process (NNGP) or thin-plate spline.
* Time: linear trend or cyclic cubic spline, optionally interacting with space.
* Responses can share fixed, random, spatial and temporal effects.

## Example

```r
library(BINTS)
sim <- simulate_bm_data(n_year = 4, n_locality = 20, n_rep = 2, seed = 1)
fit <- fit_model(sim, fixed_effects = c("x1", "x2"),
                 temporal_effect = "linear", spatial_effect = "coordinates",
                 n_chains = 2, n_iter = 2000, n_burnin = 1000, n_thin = 2,
                 parallel = FALSE)
summary(fit)
posterior_predictive_check(fit)
predict(fit)
```

For real count data use
`family = c(biomass = "truncnorm", richness = "negbinomial")`.
See `vignette("BINTS-workflow")` for a full workflow.
