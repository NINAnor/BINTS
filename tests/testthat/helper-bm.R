bm_test_data <- function(seed = 4711, ...) {
  simulate_bm_data(n_year = 5, n_locality = 25, n_rep = 4, seed = seed, ...)
}

# Short run used in automated tests
bm_test_fit <- function(data, ..., n_iter = 1500, n_burnin = 750, n_thin = 3,
                        n_chains = 2, seed = 1, parallel = FALSE) {
  suppressMessages(suppressWarnings(
    fit_model(data, ..., n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
              n_chains = n_chains, seed = seed, parallel = parallel)
  ))
}

# Fitted models are expensive (NIMBLE compiles C++); fit each once per test run
.bm_cache <- new.env()
bm_cached <- function(key, expr) {
  if (!exists(key, envir = .bm_cache, inherits = FALSE)) {
    assign(key, expr, envir = .bm_cache)
  }
  get(key, envir = .bm_cache)
}

# Fitted models shared by several test files
fit_joint <- function() bm_cached("joint", {
  d <- bm_test_data()
  list(d = d, fit = bm_test_fit(d, fixed_effects = c("x1", "x2"), n_iter = 2500,
                                n_burnin = 1250, n_thin = 3))
})
fit_spatial <- function() bm_cached("spatial", {
  d <- simulate_bm_data(n_year = 5, n_locality = 30, n_rep = 3, seed = 21,
                        spatial = "coordinates", spatial_sd = c(biomass = 6, richness = 60),
                        trend = c(biomass = 1.5, richness = 0))
  list(d = d, fit = bm_test_fit(d, responses = "biomass", fixed_effects = "x1",
                                temporal_effect = "linear", spatial_effect = "coordinates",
                                random_effects = c("year", "locality")))
})
fit_cyclic <- function() bm_cached("cyclic", {
  d <- bm_test_data(seed = 7, cyclic_amplitude = c(biomass = 6, richness = 0))
  list(d = d, fit = bm_test_fit(d, responses = "biomass", temporal_effect = "cyclic",
                                random_effects = c("year", "locality")))
})
fit_nngp <- function() bm_cached("nngp", {
  d <- simulate_bm_data(n_year = 3, n_locality = 25, n_rep = 3, seed = 22,
                        spatial = "autocorrelation", spatial_sd = c(biomass = 6, richness = 60),
                        spatial_range_km = 250)
  list(d = d, fit = bm_test_fit(d, responses = "biomass", spatial_effect = "autocorrelation",
                                random_effects = "year"))
})



# Data with negative-binomial richness: richness parameters are on the log scale
bm_nb_data <- function(seed = 5, ...) {
  simulate_bm_data(
    n_year = 4, n_locality = 30, n_rep = 3, seed = seed,
    family = c(biomass = "truncnorm", richness = "negbinomial"),
    beta0 = c(biomass = 20, richness = log(500)),
    beta = list(biomass = c(x1 = 5, x2 = 0), richness = c(x1 = 0, x2 = 0.3)),
    sd_year = c(biomass = 3, richness = 0.1), sd_locality = c(biomass = 6, richness = 0.3),
    sd_replication = c(biomass = 4, richness = 0.1), nb_size = 6, ...
  )
}
fit_nb <- function() bm_cached("nb", {
  d <- bm_nb_data()
  list(d = d, fit = bm_test_fit(d, responses = "richness", family = "negbinomial",
                                fixed_effects = "x2", n_iter = 3000, n_burnin = 1500))
})
fit_pois <- function() bm_cached("pois", {
  d <- bm_nb_data()
  list(d = d, fit = bm_test_fit(d, responses = "richness", family = "poisson",
                                fixed_effects = "x2", n_iter = 3000, n_burnin = 1500))
})
fit_counts <- function() bm_cached("counts", {
  d <- bm_nb_data()
  list(d = d, fit = bm_test_fit(d, responses = c("biomass", "richness"),
                                family = c(biomass = "lognormal", richness = "negbinomial"),
                                fixed_effects = list(biomass = "x1", richness = "x2"),
                                n_iter = 3000, n_burnin = 1500))
})
