# Plan: BINTS NNGP model, allInsectsData fit, slim vignette

## Findings
- Package is BINTS (NIMBLE). `fit_model()` in R/bm_fit.R already supports: truncnorm/lognormal biomass; poisson/negbinomial richness (R/bm_families.R); NNGP spatial field (R/bm_spatial.R, bm_distributions.R `bm_nngp_code`, m=10, exponential, phi grid of 20); cyclic spline (R/bm_temporal.R); shared effects/loadings (R/bm_sharing.R); predict, PPC, residuals, plots.
- No README content, no NEWS. `allInsectsData.R` not found in repo (path to be supplied by user).
- Decision: keep NIMBLE, tighten NNGP toward spOccupancy conventions (no sampler package fits joint count+biomass).

## Step 1: Model audit and tightening (few round trips)
1. Check NNGP vs spOccupancy: default m=15, neighbour search by nearest distance among earlier-ordered sites, phi with uniform/discrete prior, spatial variance prior inverse-gamma, latent w sampled via block updates. Change defaults/code only where they differ and matter.
2. Confirm family validation: biomass in {truncnorm, lognormal}; richness in {poisson, negbinomial}; make negbinomial/poisson the richness default instead of truncnorm.
3. Add/extend tests: families, NNGP coefficients vs exact GP, short joint fit with NB + truncnorm + NNGP + cyclic + shared. Run with devtools::test() once.

## Step 2: allInsectsData.R fit
1. Read the script at the user-given path; reuse its data prep.
2. Fit `fit_model(responses = c("biomass","richness"), family = list(biomass="truncnorm", richness="negbinomial"), spatial_effect="autocorrelation", temporal_effect="cyclic", shared_spatial_effect=TRUE, shared_temporal_effect=TRUE, shared fixed/random as appropriate, spatiotemporal_interaction as optional)`; moderate MCMC, saved to scratch/outputs (.rds).
3. Convergence diagnostics, `posterior_predictive_check()` for both responses.
4. `predict()` on a grid over the study region (sf polygon), `plot_spatial_predictions()`.
5. `residuals()` + `plot_residuals()`, `plot_observed_predicted()`; save figures to a figures folder.

## Step 3: Vignette, README, NEWS
1. Shrink vignettes/BINTS-workflow.Rmd: ~20 localities, 3 years, 2 reps, 1-2 chains, ~1500 iter, family NB + truncnorm, NNGP small m; target a few minutes; keep eval guard.
2. Write README.md (purpose, install, quick example, model overview, links) and NEWS.md (0.2.0 / next version entries).
3. Render vignette once to check runtime.

## Open item
- Path to `allInsectsData.R` still needed before Step 2.
