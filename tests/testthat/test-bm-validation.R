test_that("simulator is reproducible and respects the lower bound", {
  a <- bm_test_data(1)
  b <- bm_test_data(1)
  expect_identical(a, b)
  expect_true(all(a$sum_wet_weight >= 0))
  expect_true(all(a$no_species >= 0))
  expect_true(all(c("beta0", "beta", "sd_year", "sd_locality",
                    "sd_replication", "sigma", "trend") %in%
                    names(attr(a, "truth"))))
})

test_that("fit_model validates its arguments", {
  d <- bm_test_data()
  expect_error(fit_model(1:3), "data frame")
  expect_error(fit_model(d, responses = "foo"), "responses")
  expect_error(fit_model(d, random_effects = "site"), "random_effects")
  expect_error(fit_model(d, temporal_effect = "quadratic"), "temporal_effect")
  expect_error(fit_model(d, spatial_effect = "gp"), "spatial_effect")
  expect_error(fit_model(d, spatiotemporal_interaction = TRUE), "both a temporal")
  expect_error(fit_model(d, temporal_options = list(k = 2), temporal_effect = "cyclic"),
               "k")
  expect_error(fit_model(d, temporal_options = list(foo = 1)), "temporal_options")
  expect_error(fit_model(d, shared_fixed_effects = TRUE), "no fixed effect")
  expect_error(fit_model(d, n_chains = 0), "n_chains")
  expect_error(fit_model(d, n_iter = 100, n_burnin = 200), "n_burnin")
  expect_error(fit_model(d, n_chains = 2, chain_seeds = 1L), "chain_seeds")
  expect_error(fit_model(d, n_chains = 2, chain_seeds = c(1L, 1L)), "chain_seeds")
  expect_error(fit_model(d, family = "poisson"), "not available")
  expect_error(fit_model(d, family = "lognormal", responses = "richness"),
               "not available")
  expect_error(fit_model(d, responses = "biomass", fixed_effects = "nope"),
               "not found")
  expect_error(fit_model(d, responses = "biomass", priors = list(foo = list())),
               "priors")
  expect_error(fit_model(d, responses = "biomass",
                         priors = list(biomass = list(bad = 1))),
               "Unknown prior")
})

test_that("log-normal family rejects zeros with a clear message", {
  d <- bm_test_data()
  d$sum_wet_weight[1] <- 0
  expect_error(fit_model(d, responses = "biomass", family = "lognormal"),
               "zero value")
})

test_that("custom column mapping is honoured", {
  d <- bm_test_data()
  names(d)[names(d) == "sum_wet_weight"] <- "bm"
  cols <- bm_default_columns()
  cols$biomass <- "bm"
  prep <- bm_prepare_data(d, "biomass", list(biomass = character()),
                          c("year", "locality", "replication"),
                          c(biomass = "truncnorm"), cols)
  expect_equal(prep$y$biomass, d$bm)
  expect_equal(length(prep$index$replication$levels), 4L)
})

test_that("missing covariate values drop rows with a message", {
  d <- bm_test_data()
  d$x1[1:3] <- NA
  expect_message(
    prep <- bm_prepare_data(d, "biomass", list(biomass = "x1"),
                            "year", c(biomass = "truncnorm"),
                            bm_default_columns()),
    "3 row"
  )
  expect_equal(prep$n, nrow(d) - 3L)
})

test_that("model code puts the truncation in the likelihood", {
  d <- bm_test_data()
  fam <- c(biomass = "truncnorm")
  spec <- list(responses = "biomass", family = fam,
               fixed_effects = list(biomass = "x1"),
               random_effects = c("year", "locality"), temporal_effect = "none")
  prep <- bm_prepare_data(d, "biomass", spec$fixed_effects, spec$random_effects,
                          fam, bm_default_columns())
  built <- bm_build_model(prep, spec, bm_default_priors(prep, fam))
  txt <- gsub("\\s+", " ", paste(deparse(built$code), collapse = " "))
  expect_match(txt, "dtnorm_lower\\(mu_biomass\\[i\\], sd = sigma_biomass, lower = 0\\)")
  expect_false("u_replication_biomass" %in% built$monitors)
})
