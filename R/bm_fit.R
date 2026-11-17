#' Fit a Bayesian hierarchical model for biomass and species richness
#'
#' Builds a NIMBLE model from the requested components, compiles it, and runs
#' several independent MCMC chains. This is the entry point of the `bm_` model
#' framework; it is independent of the older `create_insect_*` functions.
#'
#' @section Model:
#' For response \eqn{r} (biomass or richness) and observation \eqn{i}, the
#' linear predictor is
#' \deqn{\mu_{ri} = \beta_{0r} + x_i^\top \beta_r + u^{year}_{r,t(i)} +
#'   u^{locality}_{r,l(i)} + u^{replication}_{r,k(i)}}
#' (plus the temporal and spatial terms described below), and the response has
#' one of the following distributions, chosen per response with `family`:
#'
#' | `family` | responses | distribution | link and second parameter |
#' |:---------|:----------|:-------------|:--------------------------|
#' | `"truncnorm"` (default) | biomass, richness | \eqn{y \sim \mathrm{N}(\mu, \sigma^2)} truncated to \eqn{[0, \infty)} | identity; `sigma_<response>` |
#' | `"lognormal"` | biomass | \eqn{\log y \sim \mathrm{N}(\mu, \sigma^2)} | `sigma_<response>`; needs `y > 0` |
#' | `"poisson"` | richness | \eqn{y \sim \mathrm{Poisson}(e^{\mu})} | log; none |
#' | `"negbinomial"` | richness | \eqn{y \sim \mathrm{NB}(\text{mean} = e^{\mu}, \text{size} = \phi)} | log; `phi_<response>` |
#'
#' For example `family = c(biomass = "truncnorm", richness = "negbinomial")`.
#' The truncation at zero of the truncated normal is part of the likelihood
#' (a custom distribution, `dtnorm_lower`, evaluated in log space), not applied
#' to predictions afterwards. The count families need non-negative integers.
#' Priors are set on the scale of the linear predictor (the response for the
#' truncated normal, its logarithm for the others; `log(1 + y)` for counts).
#'
#' @section Interpreting parameters by family:
#' *Truncated normal*: `mu` and `sigma` are the location and scale of the
#' *untruncated* normal; they are **not** the mean and sd of the response. The
#' expected response is
#' \eqn{\mu + \sigma \phi(\mu/\sigma)/\Phi(\mu/\sigma)}, which is at least `mu`
#' and exceeds it when `mu` is small relative to `sigma`. `predict()` returns this
#' truncated mean by default. *Log-normal*: `mu` is the mean of `log(y)` and
#' `exp(mu)` the median; the expected value is `exp(mu + sigma^2/2)`.
#' *Poisson and negative binomial*: `mu` is the log of the mean count, so a
#' coefficient is a log rate ratio per standard deviation of the covariate and
#' `exp(mu)` is the expected count; `type = "location"` in `predict()` returns
#' the log scale.
#'
#' @section Richness distribution:
#' Richness is a count. The truncated normal is continuous, so it ignores
#' discreteness and puts no mass on exact zeros; it is adequate when counts are
#' large but can be poor for small counts. `"poisson"` assumes the variance
#' equals the mean, which richness data rarely satisfy (the posterior predictive
#' check shows this as an extreme standard deviation and deviance).
#' `"negbinomial"` adds a size parameter `phi` (variance
#' \eqn{m + m^2/\phi}, so small `phi` means strong overdispersion; its prior is
#' `Gamma(disp_shape = 1, disp_rate = 0.1)`, mean 10). Models can be compared
#' with [posterior_predictive_check()]. The `family` registry makes further
#' families easy to add.
#'
#' @section Random effects:
#' Each random effect is `Normal(0, sd)` with a half-t(3) prior on `sd`. The
#' replication effect uses the trailing number of the replication column (see
#' [bm_default_columns()]). Random effects with few levels (such as six years)
#' give weakly identified sds.
#'
#' @section Covariates:
#' Fixed-effect covariates are standardised internally; coefficients are per
#' standard deviation. Factors and characters are dummy coded. Rows with
#' missing values in any model variable are dropped with a message.
#'
#' @section Seeds and reproducibility:
#' A master `seed` (random if `NULL`, and recorded) generates one unique seed
#' per chain, used for both initial values and MCMC. Chain seeds can instead be
#' supplied via `chain_seeds`. The seeds actually used are stored in
#' `fit$seeds`. Results are identical for the same data, specification and
#' seeds whether chains run in parallel or not.
#'
#' @section Temporal effects:
#' Structured temporal effects are added to the linear predictor in addition
#' to the random effects.
#'
#' * `"linear"`: \eqn{\beta^{trend}_r (year - \overline{year})}. The slope is
#'   per year; the year is centred so the trend is orthogonal to the
#'   intercept.
#' * `"cyclic"`: a penalised cubic cyclic regression spline \eqn{f_r(s)} of
#'   a seasonal variable \eqn{s} (default column `season`, day of year),
#'   with period `temporal_options$period` (default 365) and basis size
#'   `temporal_options$k` (default 8). The spline has a sum-to-zero constraint
#'   and is fitted as a random effect with one half-t smoothing sd.
#'
#' *Structured trend versus interannual variability.* The year random effect
#' captures unstructured year-to-year deviations (shocks, weather) with
#' mean zero; the linear trend captures a monotone change over time. When both
#' are present, the year effects are deviations around the trend line, so the
#' slope is identified by the pattern across years relative to the year sd. With
#' few years (say < 5), the year effect can absorb much of any trend and the
#' slope will be uncertain. The cyclic effect describes the *within-year*
#' (seasonal) cycle and therefore does not compete with the year effect.
#' However, if replicates correspond to sampling occasions, the replication
#' random effect can also absorb seasonal variation and is only partly
#' separable from the cyclic effect; consider dropping `"replication"` from
#' `random_effects` when fitting a cyclic effect.
#'
#' @section Spatial effects:
#' Spatial effects are defined at *sites*: localities, each with one coordinate
#' pair that is constant over years and replicates (checked). Coordinates are
#' taken from `columns$longitude`/`columns$latitude` or from `coordinates =
#' c(x, y)`. By default they are longitude/latitude in degrees
#' (`spatial_options$crs = 4326`, great-circle distances); for projected
#' coordinates give `spatial_options$crs` (assumed metres; `km_per_unit`
#' converts other units).
#'
#' * `"autocorrelation"`: a nearest-neighbour Gaussian process (NNGP) with
#'   exponential correlation `exp(-d / ell)` (practical range about `3 ell`) and
#'   marginal sd `sd_sp`. Sites are ordered by longitude and each conditions on
#'   its `m` (default 10) nearest earlier sites, so the cost grows linearly in
#'   the number of sites. The length scale `ell` is restricted to a log-spaced
#'   grid (`n_range` = 20 values, by default from the median nearest-neighbour
#'   distance to half the extent of the sites, or `spatial_options$range`
#'   in km) with a uniform prior over grid points; this avoids matrix algebra
#'   inside the sampler. Predictions at unsampled locations condition on the
#'   `m` nearest fitted sites. *Assumptions*: stationary, isotropic, exponential
#'   correlation; Gaussian field on the model (pre-truncation) scale. The
#'   current implementation indexes neighbours dynamically, whose cost grows
#'   roughly quadratically with the number of sites; it is comfortable up to
#'   a few hundred to ~1000 sites.
#' * `"coordinates"`: a thin-plate regression spline of the planar (km)
#'   coordinates, with `spatial_options$k` (default 25) basis functions. The
#'   unpenalised plane (a linear trend in the coordinates) gets vague priors;
#'   the wiggly part is a random effect with one smoothing sd. *Assumptions*:
#'   smooth surface, isotropic after conversion to km, no spatial correlation
#'   beyond the spline; extrapolation beyond the sampled area reverts towards
#'   the plane and should not be trusted. It is deterministic given its
#'   coefficients, so it is cheap and its predictions are straightforward.
#'
#' The locality random effect (an independent site-level term, i.e. a
#' *nugget*) can be kept alongside a spatial effect. The two are separated only
#' by the assumed spatial structure, so with one value per site they are
#' weakly identified; drop `"locality"` from `random_effects` for a
#' simpler, better-mixing model. The NNGP field has mean zero and the spline
#' is constrained to sum to zero over the sites, so neither competes with the
#' intercept except through the usual weak confounding with very long ranges.
#'
#' @section Spatio-temporal models:
#' Any temporal effect can be combined with any spatial effect. By default the
#' combination is *additive*: the temporal term and the spatial field are
#' simply summed in the linear predictor, so the spatial pattern is the same in
#' every year (and the temporal pattern is the same everywhere).
#'
#' With `spatiotemporal_interaction = TRUE` a second spatial field (same type
#' as `spatial_effect`, with its own sd and, for the NNGP, its own length
#' scale; parameters `*_spt_*`) modulates the temporal term, so that the
#' spatial pattern changes over time:
#'
#' * `temporal_effect = "linear"`: a *spatially varying slope*,
#'   \eqn{(\beta^{trend} + g_l)(year - \overline{year})}, where \eqn{g_l} is a
#'   mean-zero field (units: change per year). The global trend is the
#'   average slope and \eqn{g_l} the local departure from it.
#' * `temporal_effect = "cyclic"`: a *spatially varying seasonal amplitude*,
#'   \eqn{g_l f(s)}, added to \eqn{f(s)}, i.e. the seasonal curve is scaled by
#'   \eqn{1 + g_l} at site \eqn{l} (dimensionless; its shape is the same
#'   everywhere).
#'
#' The interaction is deliberately low-dimensional (one extra field rather
#' than an independent field per year). It does not represent year-specific
#' spatial patterns unrelated to the trend or seasonal curve; those would need
#' an additional model component. Identifiability: the interaction field has
#' mean zero (spline: sum-to-zero over sites), so the global trend/seasonal
#' effect keeps its meaning; for the NNGP with a linear trend a shift sampler
#' improves mixing between the global slope and the field. The interaction
#' needs both a temporal and a spatial effect and enough years to estimate
#' local slopes (at least 4-5), and is weakly informed with few observations
#' per site.
#'
#' @section Joint models: sharing parameters, processes and correlations:
#' With `responses = c("biomass", "richness")` the two responses are
#' modelled jointly. Conditional on the linear predictors the two are
#' independent, and each has its own intercept and residual sd. Everything
#' else can be specified effect by effect, and three different ideas must be
#' distinguished:
#'
#' 1. *Sharing a parameter* (`shared_loadings = FALSE`): the second response
#'    uses exactly the same coefficients or random-effect values as the first.
#'    This only makes sense when the responses are on comparable scales; the
#'    covariates are standardised but biomass and richness are not, so it
#'    is mainly useful for effects that are meant to be identical in
#'    absolute units.
#' 2. *Sharing a latent ecological process* (`shared_loadings = TRUE`, the
#'    default): a single latent effect \eqn{z} (a coefficient, a set of random
#'    effects, a trend, a seasonal curve, a spatial field) enters the first
#'    response as \eqn{z} and the second as \eqn{\lambda z}. The loading
#'    \eqn{\lambda} (`lambda_*_<second response>`, prior `Normal(0, loading_sd)`
#'    with `loading_sd` = 5 times the response-scale ratio) lets the responses
#'    react with different magnitude, or in opposite directions, to the same
#'    process. *Identifiability*: the first response in `responses` is the
#'    reference with loading fixed at 1, which fixes the sign and the scale
#'    of \eqn{z}; it must therefore respond to the effect (if its effect is
#'    near zero, \eqn{\lambda} is poorly determined and the posterior can be
#'    bimodal). Each shared effect has its own loading. The second response's
#'    random-effect sd is `abs(lambda) * sd_reference`; parameters of the
#'    second response (`beta_*`, `u_*`, `w_*`, ...) are returned as derived
#'    quantities so that prediction code is the same as for separate effects.
#' 3. *Correlated response-specific effects* (`correlated_random_effects`):
#'    each response has its own random effect with its own sd, drawn
#'    jointly from a bivariate normal with correlation `rho_*` (a uniform
#'    prior on (-1, 1)). This is more flexible than a shared process (the
#'    effects need not be proportional) but needs more data. Correlated
#'    spatial or temporal effects are not implemented; a shared field with a
#'    loading is the one-factor version of such a model.
#'
#' The shared spatial effect includes the space-time interaction field.
#' The shared temporal effect shares the trend or the seasonal curve. Sharing
#' the interaction field while sharing the trend as well mixes more slowly
#' (the shift sampler for the slope is switched off).
#'
#' @param data A data frame, one row per sampling event.
#' @param responses `"biomass"`, `"richness"`, or both.
#' @param fixed_effects `NULL`, a character vector of covariate columns used for
#'   every response, or a list named by response.
#' @param random_effects Subset of `c("year", "locality", "replication")`, or
#'   `NULL`.
#' @param temporal_effect `"none"`, `"linear"` or `"cyclic"`; see Temporal
#'   effects.
#' @param spatial_effect `"none"`, `"autocorrelation"` or `"coordinates"`; see
#'   Spatial effects.
#' @param temporal_options List with `k` (cyclic basis dimension, default 8) and
#'   `period` (default 365).
#' @param spatial_options List with `m` (neighbours, 10), `n_range` (length-scale
#'   grid size, 20), `range` (`c(min, max)` length scales in km), `k` (spline
#'   basis size, 25), `crs` (4326) and `km_per_unit` (0.001, for projected
#'   coordinates).
#' @param spatiotemporal_interaction Logical; add a spatial field that modulates
#'   the temporal term (spatially varying slope or seasonal amplitude). Requires
#'   both a temporal and a spatial effect; see Spatio-temporal models.
#' @param shared_fixed_effects `FALSE`, `TRUE` (all covariates common to both
#'   responses) or a character vector of covariate names to share.
#' @param shared_random_effects `FALSE`, `TRUE` (all) or a character vector
#'   with a subset of `random_effects` to share.
#' @param shared_spatial_effect,shared_temporal_effect `TRUE` to share the
#'   spatial (including its interaction field) or temporal latent effect.
#' @param correlated_random_effects `FALSE`, `TRUE` or a character vector of
#'   random effects for which the two responses have separate effects with a
#'   correlation `rho` (not combinable with sharing the same effect).
#' @param shared_loadings `TRUE` (default): shared effects are a shared latent
#'   process scaled by a loading `lambda` for the second response. `FALSE`:
#'   the second response uses the identical parameter values.
#' @param coordinates `NULL` (use `columns`) or the names of the x and y
#'   coordinate columns, `c("longitude", "latitude")`.
#' @param priors `NULL` or a list named by response, each a named list overriding
#'   any of `beta0_mean`, `beta0_sd`, `beta_sd`, `sd_re_scale`, `sigma_scale`,
#'   and for the negative binomial `disp_shape`, `disp_rate`.
#'   Defaults are weakly informative and relative to the response scale.
#' @param family `NULL` (truncated normal for all), a single string, or a named
#'   vector by response. Options: `"truncnorm"` (biomass or richness),
#'   `"lognormal"` (biomass only; strictly positive values), and `"poisson"` or
#'   `"negbinomial"` (richness only; non-negative integers).
#' @param columns Named list mapping model roles to data columns; see
#'   [bm_default_columns()].
#' @param n_chains,n_iter,n_burnin,n_thin MCMC settings.
#' @param seed Master seed, or `NULL` for a random (recorded) one.
#' @param chain_seeds Optional vector of `n_chains` unique integer seeds.
#' @param parallel Run chains on separate worker processes (PSOCK).
#'
#' @return An object of class `bm_fit`, a list with elements `samples` (a
#'   `coda::mcmc.list`), `diagnostics` (posterior summary with R-hat and ESS
#'   per parameter), `seeds`, `mcmc`, `priors`, `spec`, `model` (NIMBLE code,
#'   constants, data, monitors), `prep` (prepared data and indices) and `call`.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 10, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", fixed_effects = "x1",
#'                  n_chains = 2, n_iter = 2000, n_burnin = 1000, n_thin = 2,
#'                  seed = 1, parallel = FALSE)
#' fit
#' summary(fit)
#' }
#' @md
fit_model <- function(data,
                      responses = c("biomass", "richness"),
                      fixed_effects = NULL,
                      random_effects = c("year", "locality", "replication"),
                      temporal_effect = "none",
                      spatial_effect = "none",
                      temporal_options = NULL,
                      spatial_options = NULL,
                      spatiotemporal_interaction = FALSE,
                      shared_fixed_effects = FALSE,
                      shared_random_effects = FALSE,
                      shared_spatial_effect = FALSE,
                      shared_temporal_effect = FALSE,
                      correlated_random_effects = FALSE,
                      shared_loadings = TRUE,
                      coordinates = NULL,
                      priors = NULL,
                      family = NULL,
                      columns = bm_default_columns(),
                      n_chains = 4,
                      n_iter = 50000,
                      n_burnin = 25000,
                      n_thin = 10,
                      seed = NULL,
                      chain_seeds = NULL,
                      parallel = TRUE) {
  cl <- match.call()
  if (!is.null(coordinates)) {
    if (!is.character(coordinates) || length(coordinates) != 2L) {
      stop("`coordinates` must be two column names: c(x, y), e.g. ",
           "c(\"longitude\", \"latitude\").", call. = FALSE)
    }
    columns$longitude <- coordinates[1]
    columns$latitude <- coordinates[2]
  }
  shared <- list(
    shared_fixed_effects = shared_fixed_effects,
    shared_random_effects = shared_random_effects,
    shared_spatial_effect = shared_spatial_effect,
    shared_temporal_effect = shared_temporal_effect,
    correlated_random_effects = correlated_random_effects
  )
  bm_validate_args(data, responses, fixed_effects, random_effects,
                   temporal_effect, spatial_effect, spatiotemporal_interaction,
                   shared, coordinates, columns, n_chains, n_iter, n_burnin,
                   n_thin, seed, chain_seeds, parallel)
  fixed <- bm_resolve_fixed(fixed_effects, responses)
  family <- bm_resolve_family(family, responses)
  sharing <- bm_resolve_sharing(
    responses, fixed, random_effects, temporal_effect, spatial_effect,
    shared_fixed_effects, shared_random_effects, shared_spatial_effect,
    shared_temporal_effect, correlated_random_effects, shared_loadings
  )

  prep <- bm_prepare_data(data, responses, fixed, random_effects, family, columns,
                          temporal_effect, temporal_options, spatial_effect,
                          spatial_options)
  if (temporal_effect == "cyclic" && "replication" %in% random_effects) {
    message("Note: a cyclic seasonal effect and a replication random effect ",
            "can be partly confounded when replicates correspond to sampling ",
            "occasions; see ?fit_model.")
  }
  priors <- bm_resolve_priors(priors, bm_default_priors(prep, family))
  spec <- list(
    responses = responses, family = family, fixed_effects = fixed,
    random_effects = random_effects, temporal_effect = temporal_effect,
    spatial_effect = spatial_effect, temporal_options = temporal_options,
    spatial_options = spatial_options,
    spatiotemporal_interaction = spatiotemporal_interaction,
    shared = shared, sharing = sharing, columns = columns
  )
  built <- bm_build_model(prep, spec, priors)

  seeds <- bm_chain_seeds(seed, chain_seeds, n_chains)
  inits <- lapply(seeds$chains, function(s) bm_make_inits(prep, spec, priors, s))

  t0 <- Sys.time()
  samples <- bm_run_mcmc(built, inits, seeds$chains, n_iter, n_burnin, n_thin,
                         parallel)
  runtime <- difftime(Sys.time(), t0, units = "secs")

  diagnostics <- bm_diagnostics(samples)
  bm_convergence_warnings(diagnostics)

  structure(
    list(
      samples = samples, diagnostics = diagnostics,
      seeds = list(master = seeds$master, chains = seeds$chains),
      mcmc = list(n_chains = n_chains, n_iter = n_iter, n_burnin = n_burnin,
                  n_thin = n_thin, parallel = parallel, runtime = runtime),
      priors = priors, spec = spec, model = built, prep = prep, call = cl
    ),
    class = "bm_fit"
  )
}
