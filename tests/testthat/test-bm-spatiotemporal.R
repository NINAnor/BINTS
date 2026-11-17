st_data <- function(spatial, n_locality = 25, n_year = 3, ...) {
  simulate_bm_data(
    n_year = n_year, n_locality = n_locality, n_rep = 3, seed = 91, spatial = spatial,
    spatial_sd = c(biomass = 6, richness = 60), spatial_range_km = 250,
    trend = c(biomass = 1.5, richness = 0),
    cyclic_amplitude = c(biomass = 0, richness = 0), ...
  )
}

build_st <- function(d, te, se, interaction) {
  fam <- c(biomass = "truncnorm")
  spec <- list(responses = "biomass", family = fam,
               fixed_effects = list(biomass = character()),
               random_effects = "year", temporal_effect = te, spatial_effect = se,
               spatiotemporal_interaction = interaction)
  prep <- bm_prepare_data(d, "biomass", spec$fixed_effects, "year", fam,
                          bm_default_columns(), te, spatial_effect = se)
  bm_build_model(prep, spec, bm_default_priors(prep, fam))
}
txt <- function(b) gsub("\\s+", " ", paste(deparse(b$code), collapse = " "))

test_that("all temporal x spatial combinations build, additive and interacting", {
  d <- st_data("autocorrelation")
  for (te in c("linear", "cyclic")) {
    for (se in c("autocorrelation", "coordinates")) {
      add <- txt(build_st(d, te, se, FALSE))
      expect_match(add, if (te == "linear") "beta_trend_biomass \\* year_c" else "f_cyc_biomass\\[i\\]")
      expect_match(add, "sd_sp_biomass")
      expect_false(grepl("spt", add))
      int <- txt(build_st(d, te, se, TRUE))
      expect_match(int, "sd_spt_biomass")
      mult <- if (te == "linear") "year_c\\[i\\]" else "f_cyc_biomass\\[i\\]"
      field <- if (se == "autocorrelation") "w_spt_biomass" else "f_spt_biomass"
      expect_match(int, paste0(field, "\\[site_id\\[i\\]\\] \\* ", mult))
    }
  }
})

test_that("interaction requires both a temporal and a spatial effect", {
  d <- st_data("autocorrelation")
  expect_error(fit_model(d, spatiotemporal_interaction = TRUE,
                         temporal_effect = "linear"), "both a temporal")
  expect_error(fit_model(d, spatiotemporal_interaction = TRUE,
                         spatial_effect = "coordinates"), "both a temporal")
  expect_error(fit_model(d, spatiotemporal_interaction = "yes"), "TRUE or FALSE")
})

test_that("interaction uses separate parameters and a trend shift sampler", {
  d <- st_data("autocorrelation")
  b <- build_st(d, "linear", "autocorrelation", TRUE)
  expect_true(all(c("w_sp_biomass", "w_spt_biomass", "phi_idx_sp_biomass",
                    "phi_idx_spt_biomass") %in% b$monitors))
  expect_true(any(vapply(b$shifts, function(s) s$target == "beta_trend_biomass",
                         logical(1))))
  expect_false(any(vapply(b$blocks, function(x) "beta_trend_biomass" %in% x, logical(1))))
  b2 <- build_st(d, "cyclic", "coordinates", TRUE)
  expect_true(all(c("v_sp_biomass", "v_spt_biomass") %in% b2$monitors))
  expect_length(b2$shifts, 0L)
})

test_that("spatio-temporal models compile and run (additive and interaction)", {
  skip_if_not_installed("nimble")
  d <- st_data("autocorrelation")
  d$sum_wet_weight <- d$sum_wet_weight + 0
  grid <- expand.grid(te = c("linear", "cyclic"), se = c("autocorrelation", "coordinates"),
                      int = c(FALSE, TRUE), stringsAsFactors = FALSE)
  for (k in seq_len(nrow(grid))) {
    g <- grid[k, ]
    fit <- bm_test_fit(d, responses = "biomass", temporal_effect = g$te,
                       spatial_effect = g$se, spatiotemporal_interaction = g$int,
                       random_effects = "year", n_iter = 1000, n_burnin = 500,
                       n_thin = 2)
    expect_s3_class(fit, "bm_fit")
    pars <- fit$diagnostics$parameter
    expect_equal("sd_spt_biomass" %in% pars, g$int,
                 info = paste(g$te, g$se, g$int))
    expect_true(all(is.finite(fit$diagnostics$mean)))
  }
})

test_that("joint responses fit with temporal + spatial + interaction", {
  skip_if_not_installed("nimble")
  d <- st_data("coordinates")
  fit <- bm_test_fit(d, responses = c("biomass", "richness"),
                     temporal_effect = "linear", spatial_effect = "coordinates",
                     spatiotemporal_interaction = TRUE, random_effects = "year",
                     n_iter = 1000, n_burnin = 500, n_thin = 2)
  expect_true(all(c("sd_spt_biomass", "sd_spt_richness") %in% fit$diagnostics$parameter))
})

test_that("linear trend + NNGP (additive) recovers trend and field (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(
    n_year = 6, n_locality = 60, n_rep = 3, seed = 11, spatial = "autocorrelation",
    spatial_sd = c(biomass = 8, richness = 80), spatial_range_km = 250,
    trend = c(biomass = 1.5, richness = 0), sd_year = c(biomass = 1, richness = 10),
    sd_locality = c(biomass = 0.01, richness = 1),
    sd_replication = c(biomass = 0.01, richness = 1)
  )
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     temporal_effect = "linear", spatial_effect = "autocorrelation",
                     random_effects = "year", n_chains = 3, n_iter = 20000,
                     n_burnin = 10000, n_thin = 10)
  est <- fit$diagnostics
  row <- est[est$parameter == "beta_trend_biomass", ]
  expect_true(row$q2.5 <= 1.5 && 1.5 <= row$q97.5)
  sites <- fit$prep$spatial$site_names
  w <- colMeans(as.matrix(fit$samples)[, paste0("w_sp_biomass[", seq_along(sites), "]")])
  expect_gt(cor(w, truth$spatial_field$biomass[sites]), 0.9)
})

test_that("space-varying trend is recovered with coordinate splines (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(
    n_year = 6, n_locality = 70, n_rep = 3, seed = 12, spatial = "coordinates",
    spatial_sd = c(biomass = 6, richness = 60),
    trend = c(biomass = 1.5, richness = 0), interaction = TRUE,
    interaction_sd = c(biomass = 1.2, richness = 0.1),
    sd_year = c(biomass = 1, richness = 10),
    sd_locality = c(biomass = 0.01, richness = 1),
    sd_replication = c(biomass = 0.01, richness = 1)
  )
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     temporal_effect = "linear", spatial_effect = "coordinates",
                     spatiotemporal_interaction = TRUE, random_effects = "year",
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  sites <- fit$prep$spatial$site_names
  est <- colMeans(bm_spatial_draws(fit, "biomass", fit$prep$spatial$xy, tag = "spt"))
  expect_gt(cor(est, truth$interaction_field$biomass[sites]), 0.8)
  row <- fit$diagnostics[fit$diagnostics$parameter == "beta_trend_biomass", ]
  expect_true(row$q2.5 <= 1.5 && 1.5 <= row$q97.5)
})

test_that("space-varying trend is recovered with an NNGP (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(
    n_year = 6, n_locality = 70, n_rep = 3, seed = 13, spatial = "autocorrelation",
    spatial_sd = c(biomass = 6, richness = 60), spatial_range_km = 300,
    trend = c(biomass = 1.5, richness = 0), interaction = TRUE,
    interaction_sd = c(biomass = 1.2, richness = 0.1),
    sd_year = c(biomass = 1, richness = 10),
    sd_locality = c(biomass = 0.01, richness = 1),
    sd_replication = c(biomass = 0.01, richness = 1)
  )
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     temporal_effect = "linear", spatial_effect = "autocorrelation",
                     spatiotemporal_interaction = TRUE, random_effects = "year",
                     n_chains = 3, n_iter = 25000, n_burnin = 12500, n_thin = 10)
  sites <- fit$prep$spatial$site_names
  w <- colMeans(as.matrix(fit$samples)[, paste0("w_spt_biomass[", seq_along(sites), "]")])
  expect_gt(cor(w, truth$interaction_field$biomass[sites]), 0.7)
  row <- fit$diagnostics[fit$diagnostics$parameter == "beta_trend_biomass", ]
  expect_true(row$q2.5 <= 1.5 && 1.5 <= row$q97.5)
})
