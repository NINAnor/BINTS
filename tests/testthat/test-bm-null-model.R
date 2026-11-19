test_that("null model compiles, runs and returns a bm_fit", {
  skip_if_not_installed("nimble")
  d <- bm_test_data()
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"))
  expect_s3_class(fit, "bm_fit")
  expect_s3_class(fit$samples, "mcmc.list")
  expect_length(fit$samples, 2L)
  expect_equal(length(fit$seeds$chains), 2L)
  expect_true(all(c("beta0_biomass", "sigma_biomass", "sd_year_biomass") %in%
                    fit$diagnostics$parameter))
  expect_output(print(fit), "bm_fit")
  expect_s3_class(summary(fit), "summary.bm_fit")
  expect_false(any(grepl("^u_", summary(fit)$table$parameter)))
  expect_true(any(grepl("^u_", summary(fit, random_levels = TRUE)$table$parameter)))
})

test_that("chain seeds are recorded, reproducible and can be supplied", {
  skip_if_not_installed("nimble")
  d <- bm_test_data()
  a <- bm_test_fit(d, responses = "biomass", seed = 11)
  b <- bm_test_fit(d, responses = "biomass", seed = 11)
  c <- bm_test_fit(d, responses = "biomass", seed = 12)
  expect_identical(a$seeds, b$seeds)
  expect_equal(as.matrix(a$samples), as.matrix(b$samples))
  expect_false(isTRUE(all.equal(as.matrix(a$samples), as.matrix(c$samples))))
  expect_equal(anyDuplicated(a$seeds$chains), 0L)

  s <- bm_test_fit(d, responses = "biomass", chain_seeds = a$seeds$chains)
  expect_equal(as.matrix(s$samples), as.matrix(a$samples))
  expect_equal(s$seeds$chains, a$seeds$chains)
})

test_that("results are identical for parallel and serial execution", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  skip_if(parallel::detectCores() < 3, "needs >= 3 cores")
  d <- bm_test_data()
  ser <- bm_test_fit(d, responses = "biomass", seed = 21)
  par <- bm_test_fit(d, responses = "biomass", seed = 21, parallel = TRUE)
  expect_equal(as.matrix(par$samples), as.matrix(ser$samples))
})

test_that("joint response model without sharing fits both responses", {
  skip_if_not_installed("nimble")
  d <- bm_test_data()
  fit <- bm_test_fit(d, responses = c("biomass", "richness"),
                     fixed_effects = list(biomass = "x1", richness = "x2"))
  pars <- fit$diagnostics$parameter
  expect_true(all(c("beta0_biomass", "beta0_richness", "sigma_richness") %in% pars))
  expect_true("beta_richness[1]" %in% pars)
})

test_that("log-normal biomass runs", {
  skip_if_not_installed("nimble")
  d <- bm_test_data()
  d$sum_wet_weight <- pmax(d$sum_wet_weight, 0.5)
  fit <- bm_test_fit(d, responses = "biomass", family = "lognormal")
  expect_s3_class(fit, "bm_fit")
})

test_that("null model recovers simulated parameters (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(n_year = 6, n_locality = 60, n_rep = 5, seed = 2024)
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  est <- fit$diagnostics
  get <- function(p) est[est$parameter == p, ]
  covered <- function(p, v) get(p)$q2.5 <= v && v <= get(p)$q97.5
  expect_true(covered("beta_biomass[1]", truth$beta$biomass[["x1"]]))
  expect_true(covered("beta_biomass[2]", truth$beta$biomass[["x2"]]))
  expect_true(covered("sigma_biomass", truth$sigma[["biomass"]]))
  expect_true(covered("sd_locality_biomass", truth$sd_locality[["biomass"]]))
  expect_true(covered("beta0_biomass", truth$beta0[["biomass"]]))
})
