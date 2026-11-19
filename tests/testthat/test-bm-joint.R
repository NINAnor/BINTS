both <- c("biomass", "richness")
fam2 <- c(biomass = "truncnorm", richness = "truncnorm")

build_joint <- function(d, fixed = list(biomass = c("x1", "x2"), richness = c("x1", "x2")),
                        random = c("year", "locality"), te = "none", se = "none",
                        interaction = FALSE, share = list()) {
  sh <- bm_resolve_sharing(
    both, fixed, random, te, se,
    shared_fixed_effects = share$fixed %||% FALSE,
    shared_random_effects = share$random %||% FALSE,
    shared_spatial_effect = share$spatial %||% FALSE,
    shared_temporal_effect = share$temporal %||% FALSE,
    correlated_random_effects = share$correlated %||% FALSE,
    shared_loadings = share$loadings %||% TRUE
  )
  spec <- list(responses = both, family = fam2, fixed_effects = fixed,
               random_effects = random, temporal_effect = te, spatial_effect = se,
               spatiotemporal_interaction = interaction, sharing = sh)
  prep <- bm_prepare_data(d, both, fixed, random, fam2, bm_default_columns(), te,
                          spatial_effect = se)
  pri <- bm_default_priors(prep, fam2)
  list(built = bm_build_model(prep, spec, pri), spec = spec, prep = prep, pri = pri)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
code_txt <- function(b) gsub("\\s+", " ", paste(deparse(b$built$code), collapse = " "))

test_that("sharing structure appears in the generated NIMBLE code", {
  d <- bm_test_data()
  b <- build_joint(d, random = "year",
                   share = list(fixed = "x1", random = TRUE))
  txt <- code_txt(b)
  expect_match(txt, "lambda_beta1_richness ~ dnorm")
  expect_match(txt, "beta_richness\\[1L\\] <- lambda_beta1_richness \\* beta_biomass\\[1L\\]")
  expect_match(txt, "beta_richness\\[2L\\] ~ dnorm")
  expect_match(txt, "u_year_richness\\[j\\] <- lambda_re_year_richness \\* u_year_biomass\\[j\\]")
  expect_match(txt, "sd_year_richness <- abs\\(lambda_re_year_richness\\) \\* sd_year_biomass")
  expect_false(grepl("sd_year_richness ~", txt))
  expect_true(all(c("lambda_beta1_richness", "lambda_re_year_richness") %in%
                    b$built$monitors))
  # derived nodes get no initial values; loadings do
  inits <- bm_make_inits(b$prep, b$spec, b$pri, 1)
  expect_true(all(c("lambda_beta1_richness", "lambda_re_year_richness") %in% names(inits)))
  expect_false("u_year_richness" %in% names(inits))
  expect_false("sd_year_richness" %in% names(inits))
  expect_true(is.na(inits$beta_richness[1]) && !is.na(inits$beta_richness[2]))
})

test_that("identical sharing has no loadings; correlated effects have rho", {
  d <- bm_test_data()
  b <- build_joint(d, fixed = list(biomass = "x1", richness = "x1"), random = "year",
                   share = list(fixed = "x1", random = TRUE, loadings = FALSE))
  txt <- code_txt(b)
  expect_match(txt, "u_year_richness\\[j\\] <- u_year_biomass\\[j\\]")
  expect_match(txt, "sd_year_richness <- sd_year_biomass")
  expect_false(grepl("lambda", txt))
  b2 <- build_joint(d, fixed = list(biomass = "x1", richness = "x1"),
                    random = c("year", "locality"),
                    share = list(correlated = "locality"))
  txt2 <- code_txt(b2)
  expect_match(txt2, "rho_locality_richness ~ dunif\\(-1, 1\\)")
  expect_match(txt2, "u_locality_richness\\[j\\] ~ dnorm\\(rho_locality_richness \\* sd_locality_richness/sd_locality_biomass")
  expect_true("rho_locality_richness" %in% b2$built$monitors)
  expect_equal(bm_make_inits(b2$prep, b2$spec, b2$pri, 1)$rho_locality_richness, 0)
})

test_that("shared spatial and temporal latents are derived", {
  d <- simulate_bm_data(n_year = 3, n_locality = 20, n_rep = 3, seed = 3,
                        spatial = "autocorrelation")
  b <- build_joint(d, random = "year", te = "linear", se = "autocorrelation",
                   share = list(spatial = TRUE, temporal = TRUE))
  txt <- code_txt(b)
  expect_match(txt, "w_sp_richness\\[1:L_site\\] <- lambda_sp_richness \\* w_sp_biomass\\[1:L_site\\]")
  expect_match(txt, "beta_trend_richness <- lambda_trend_richness \\* beta_trend_biomass")
  expect_match(txt, "phi_idx_sp_richness <- phi_idx_sp_biomass")
  # reference shift sampler compensates through the second intercept
  sh1 <- b$built$shifts[[1]]
  expect_equal(sh1$target, "beta0_biomass")
  expect_equal(sh1$target2, "beta0_richness")
  expect_equal(sh1$lamNode, "lambda_sp_richness")
  # no shift (and no own field sampler) for the derived second response
  expect_false(any(grepl("w_sp_richness", unlist(lapply(b$built$shifts, `[[`, "nodes")))))
  b2 <- build_joint(d, random = "year", te = "linear", se = "autocorrelation",
                    interaction = TRUE, share = list(spatial = TRUE))
  tr <- Filter(function(s) s$target == "beta_trend_biomass", b2$built$shifts)[[1]]
  expect_equal(tr$target2, "beta_trend_richness")
  expect_equal(tr$lamNode, "lambda_spt_richness")
  # sharing both slope and interaction field: slope shift is switched off
  b3 <- build_joint(d, random = "year", te = "linear", se = "autocorrelation",
                    interaction = TRUE, share = list(spatial = TRUE, temporal = TRUE))
  expect_false(any(vapply(b3$built$shifts, function(s) s$target == "beta_trend_biomass",
                          logical(1))))
})

test_that("sharing arguments are validated", {
  d <- bm_test_data()
  f <- function(...) fit_model(d, ...)
  expect_error(f(shared_random_effects = "locality", random_effects = "year"),
               "not among")
  expect_error(f(shared_spatial_effect = TRUE), "needs a spatial")
  expect_error(f(shared_temporal_effect = TRUE), "needs a temporal")
  expect_error(f(shared_fixed_effects = "x1", fixed_effects = list(biomass = "x1")),
               "not among")
  expect_error(f(shared_fixed_effects = TRUE), "no fixed effect")
  expect_error(f(shared_random_effects = "year", correlated_random_effects = "year"),
               "both shared and correlated")
  expect_error(f(shared_loadings = NA), "shared_loadings")
  expect_error(f(shared_random_effects = 3), "character vector")
})

test_that("joint models with each kind of sharing compile and run", {
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(n_year = 4, n_locality = 25, n_rep = 3, seed = 5,
                        spatial = "autocorrelation", trend = c(biomass = 1, richness = 10),
                        loadings = c(spatial = 10, trend = 10, year = 15))
  run <- function(...) bm_test_fit(d, responses = both, n_iter = 1000, n_burnin = 500,
                                   n_thin = 2, ...)
  f1 <- run(fixed_effects = c("x1", "x2"), shared_fixed_effects = "x1",
            shared_random_effects = TRUE, random_effects = c("year", "locality"))
  expect_true(all(c("lambda_beta1_richness", "lambda_re_year_richness",
                    "lambda_re_locality_richness") %in% f1$diagnostics$parameter))
  expect_output(print(f1), "Shared")
  f2 <- run(shared_random_effects = "year", random_effects = c("year", "locality"),
            shared_loadings = FALSE)
  expect_false(any(grepl("^lambda", f2$diagnostics$parameter)))
  f3 <- run(random_effects = c("year", "locality"), correlated_random_effects = "locality")
  expect_true("rho_locality_richness" %in% f3$diagnostics$parameter)
  f4 <- run(random_effects = "year", temporal_effect = "linear",
            spatial_effect = "autocorrelation", shared_spatial_effect = TRUE,
            shared_temporal_effect = TRUE, shared_random_effects = "year")
  expect_true(all(c("lambda_sp_richness", "lambda_trend_richness") %in%
                    f4$diagnostics$parameter))
  f5 <- run(random_effects = "year", temporal_effect = "linear",
            spatial_effect = "coordinates", spatiotemporal_interaction = TRUE,
            shared_spatial_effect = TRUE)
  expect_true(all(c("lambda_sp_richness", "lambda_spt_richness") %in%
                    f5$diagnostics$parameter))
  # derived nodes are reproduced exactly from their definition
  m <- as.matrix(f4$samples)
  expect_equal(m[, "u_year_richness[2]"],
               m[, "lambda_re_year_richness"] * m[, "u_year_biomass[2]"])
  expect_equal(m[, "beta_trend_richness"],
               m[, "lambda_trend_richness"] * m[, "beta_trend_biomass"])
})

test_that("shared fixed effect and year effect loadings are recovered (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(n_year = 8, n_locality = 50, n_rep = 3, seed = 31,
                        loadings = c(x1 = 30, year = 20),
                        sd_year = c(biomass = 3, richness = 0),
                        sd_replication = c(biomass = 0.01, richness = 0.01))
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = both, fixed_effects = c("x1", "x2"),
                     random_effects = c("year", "locality"),
                     shared_fixed_effects = "x1", shared_random_effects = "year",
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  est <- fit$diagnostics
  covered <- function(p, v) {
    r <- est[est$parameter == p, ]
    r$q2.5 <= v && v <= r$q97.5
  }
  expect_true(covered("lambda_beta1_richness", 30))
  expect_true(covered("lambda_re_year_richness", 20))
  expect_true(covered("beta_biomass[1]", truth$beta$biomass[["x1"]]))
  expect_true(covered("beta_richness[1]", 30 * truth$beta$biomass[["x1"]]))
  expect_true(covered("sigma_richness", truth$sigma[["richness"]]))
})

test_that("correlation between response-specific random effects is recovered (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(n_year = 3, n_locality = 120, n_rep = 3, seed = 32,
                        rho = c(locality = 0.7),
                        sd_year = c(biomass = 0.01, richness = 0.01),
                        sd_replication = c(biomass = 0.01, richness = 0.01),
                        sd_locality = c(biomass = 6, richness = 100))
  fit <- bm_test_fit(d, responses = both, random_effects = "locality",
                     correlated_random_effects = "locality",
                     n_chains = 3, n_iter = 15000, n_burnin = 7500, n_thin = 10)
  # The model estimates the correlation realised among the simulated sites,
  # which differs from the nominal 0.7 by sampling error
  te <- attr(d, "truth")$locality_effect
  realised <- stats::cor(te$biomass, te$richness)
  r <- fit$diagnostics[fit$diagnostics$parameter == "rho_locality_richness", ]
  expect_true(r$q2.5 <= realised && realised <= r$q97.5)
  expect_gt(r$mean, 0.4)
})

test_that("a shared spatial field with loading is recovered (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(n_year = 3, n_locality = 60, n_rep = 3, seed = 33,
                        spatial = "autocorrelation", spatial_sd = c(biomass = 6, richness = 0),
                        spatial_range_km = 250, loadings = c(spatial = 12),
                        sd_locality = c(biomass = 0.01, richness = 1),
                        sd_replication = c(biomass = 0.01, richness = 1),
                        sd_year = c(biomass = 1, richness = 10))
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = both, random_effects = "year",
                     spatial_effect = "autocorrelation", shared_spatial_effect = TRUE,
                     n_chains = 3, n_iter = 25000, n_burnin = 12500, n_thin = 10)
  r <- fit$diagnostics[fit$diagnostics$parameter == "lambda_sp_richness", ]
  expect_true(r$q2.5 <= 12 && 12 <= r$q97.5)
  sites <- fit$prep$spatial$site_names
  wr <- colMeans(as.matrix(fit$samples)[, paste0("w_sp_richness[", seq_along(sites), "]")])
  expect_gt(cor(wr, truth$spatial_field$richness[sites]), 0.85)
})
