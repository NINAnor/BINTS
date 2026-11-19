test_that("temporal terms enter the model code only when requested", {
  d <- bm_test_data()
  fam <- c(biomass = "truncnorm")
  mk <- function(te) {
    spec <- list(responses = "biomass", family = fam,
                 fixed_effects = list(biomass = character()),
                 random_effects = "year", temporal_effect = te)
    prep <- bm_prepare_data(d, "biomass", spec$fixed_effects, "year", fam,
                            bm_default_columns(), te)
    bm_build_model(prep, spec, bm_default_priors(prep, fam))
  }
  txt <- function(b) gsub("\\s+", " ", paste(deparse(b$code), collapse = " "))
  expect_false(grepl("beta_trend|f_cyc", txt(mk("none"))))
  lin <- mk("linear")
  expect_match(txt(lin), "beta_trend_biomass \\* year_c\\[i\\]")
  expect_equal(mean(lin$constants$year_c), 0)
  cyc <- mk("cyclic")
  expect_match(txt(cyc), "Bcyc\\[1:n, 1:K_cyc\\] %\\*% w_cyc_biomass\\[1:K_cyc\\]")
  expect_equal(ncol(cyc$constants$Bcyc), cyc$constants$K_cyc)
})

test_that("cyclic basis is centred and validates the seasonal range", {
  cb <- bm_cyclic_basis(1:364, k = 8, period = 365)
  expect_equal(unname(colSums(cb$basis)), rep(0, cb$K), tolerance = 1e-8)
  expect_error(bm_cyclic_basis(c(1, 400), 8, 365), "period")
  d <- bm_test_data()
  d$start_date_julian <- d$start_date_julian + 300L
  expect_error(fit_model(d, responses = "biomass", temporal_effect = "cyclic"),
               "period")
})

test_that("linear temporal model compiles and runs", {
  skip_if_not_installed("nimble")
  d <- bm_test_data()
  fit <- bm_test_fit(d, responses = c("biomass", "richness"),
                     temporal_effect = "linear")
  expect_true(all(c("beta_trend_biomass", "beta_trend_richness") %in%
                    fit$diagnostics$parameter))
})

test_that("cyclic temporal model compiles, runs and gives cyclic draws", {
  skip_if_not_installed("nimble")
  d <- bm_test_data()
  fit <- bm_test_fit(d, responses = "biomass", temporal_effect = "cyclic",
                     random_effects = c("year", "locality"))
  expect_true("sd_cyc_biomass" %in% fit$diagnostics$parameter)
  dr <- bm_cyclic_draws(fit, "biomass", c(100, 200, 300))
  expect_equal(dim(dr), c(nrow(as.matrix(fit$samples)), 3L))
  # cyclic: day 0 and day 365 are the same point
  ends <- bm_cyclic_draws(fit, "biomass", c(0, 365))
  expect_equal(ends[, 1], ends[, 2], tolerance = 1e-8)
  expect_false(any(grepl("w_cyc", summary(fit)$table$parameter)))
})

test_that("linear trend is recovered (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(n_year = 8, n_locality = 40, n_rep = 4, seed = 818,
                        trend = c(biomass = 1.5, richness = -25))
  fit <- bm_test_fit(d, fixed_effects = NULL, temporal_effect = "linear",
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  est <- fit$diagnostics
  tr <- attr(d, "truth")$trend
  for (r in c("biomass", "richness")) {
    row <- est[est$parameter == paste0("beta_trend_", r), ]
    expect_true(row$q2.5 <= tr[[r]] && tr[[r]] <= row$q97.5, info = r)
  }
})

test_that("cyclic effect is recovered (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  # No replication effect in the truth: replicates are sampling occasions and
  # would otherwise be partly confounded with the seasonal effect.
  d <- simulate_bm_data(n_year = 6, n_locality = 40, n_rep = 6, seed = 303,
                        sd_replication = c(biomass = 0.01, richness = 0.01),
                        cyclic_amplitude = c(biomass = 8, richness = 0))
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     temporal_effect = "cyclic",
                     random_effects = c("year", "locality"),
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  doy <- seq(125, 265, by = 10)
  est <- colMeans(bm_cyclic_draws(fit, "biomass", doy))
  truth <- attr(d, "truth")$cyclic_fun(doy, "biomass")
  expect_gt(cor(est, truth), 0.95)
  expect_lt(sqrt(mean(((est - mean(est)) - (truth - mean(truth)))^2)), 1)
})
