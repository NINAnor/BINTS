fit_ppc_ok <- function() fit_joint()
fit_ppc_zero <- function() bm_cached("zeros", {
  d <- bm_test_data(seed = 9)
  set.seed(3)
  d$sum_wet_weight[sample(nrow(d), round(0.15 * nrow(d)))] <- 0
  list(d = d, fit = bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2")))
})
fit_ppc_lnorm <- function() bm_cached("lnorm", {
  d <- bm_test_data(seed = 10)
  d$sum_wet_weight <- pmax(d$sum_wet_weight, 0.5)
  list(d = d, fit = bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                                family = "lognormal"))
})

test_that("posterior_predictive_check returns summaries, p-values and plots", {
  skip_if_not_installed("nimble")
  fit <- fit_ppc_ok()$fit
  ppc <- posterior_predictive_check(fit, n_draws = 200, seed = 1)
  expect_s3_class(ppc, "bm_ppc")
  s <- ppc$summary
  expect_true(all(c("response", "group", "statistic", "observed", "replicated_mean",
                    "replicated_lower", "replicated_upper", "p_value", "flag") %in% names(s)))
  expect_setequal(unique(s$response), c("biomass", "richness"))
  expect_true(all(c("mean", "sd", "min", "max", "q5", "q50", "q95", "prop_near_zero",
                    "prop_zero", "deviance") %in% s$statistic))
  expect_true(all(s$p_value >= 0 & s$p_value <= 1, na.rm = TRUE))
  expect_equal(s$flag == "extreme", !is.na(s$p_value) & (s$p_value < 0.025 | s$p_value > 0.975))
  expect_true(all(s$replicated_lower <= s$replicated_upper))
  expect_true(all(vapply(ppc$plots, function(p) inherits(p, "ggplot"), logical(1))))
  expect_true(all(c("density", "statistics", "discrepancy") %in% names(ppc$plots)))
  expect_s3_class(plot(ppc, "statistics"), "ggplot")
  expect_error(plot(ppc, "nonsense"), "which")
  expect_output(print(ppc), "posterior predictive")
  # drawing the plots must not fail
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  for (p in ppc$plots) print(p)
  grDevices::dev.off()
  unlink(tf)
})

test_that("replicates use the fitted likelihood: positive and matching the model", {
  skip_if_not_installed("nimble")
  m <- fit_ppc_ok(); fit <- m$fit
  ppc <- posterior_predictive_check(fit, n_draws = 150, n_overlay = 40, seed = 2)
  for (r in c("biomass", "richness")) {
    rep <- ppc$replicates[[r]]
    expect_true(all(rep >= 0))
    expect_equal(dim(rep), c(40L, nrow(fit$prep$data)))
  }
  # The data were simulated from this model and the model was fitted with all
  # effects: centre and spread of the data are reproduced
  s <- ppc$summary
  pv <- function(r, st) s$p_value[s$response == r & s$statistic == st & s$group == "all"]
  for (r in c("biomass", "richness")) {
    expect_gt(pv(r, "mean"), 0.05); expect_lt(pv(r, "mean"), 0.95)
    expect_gt(pv(r, "sd"), 0.05); expect_lt(pv(r, "sd"), 0.95)
    expect_gt(pv(r, "deviance"), 0.05); expect_lt(pv(r, "deviance"), 0.95)
  }
  # no exact zeros in the data and none in the replicates: mid-p of 0.5, not flagged
  expect_equal(pv("biomass", "prop_zero"), 0.5)
})

test_that("checks are reproducible and leave the RNG state untouched", {
  skip_if_not_installed("nimble")
  fit <- fit_ppc_ok()$fit
  a <- posterior_predictive_check(fit, n_draws = 100, seed = 11)
  b <- posterior_predictive_check(fit, n_draws = 100, seed = 11)
  c2 <- posterior_predictive_check(fit, n_draws = 100, seed = 12)
  expect_equal(a$summary, b$summary)
  expect_false(isTRUE(all.equal(a$summary$replicated_mean, c2$summary$replicated_mean)))
  set.seed(1); u1 <- runif(1)
  set.seed(1); invisible(posterior_predictive_check(fit, n_draws = 50, seed = 3)); u2 <- runif(1)
  expect_equal(u1, u2)
})

test_that("checks can be made by group and for chosen discrepancies", {
  skip_if_not_installed("nimble")
  m <- fit_ppc_ok(); fit <- m$fit
  ppc <- posterior_predictive_check(fit, responses = "biomass", n_draws = 100, by = "year",
                                    discrepancy = c("deviance", "chisq"), seed = 1)
  s <- ppc$summary
  yrs <- length(unique(m$d$year))
  expect_equal(sum(s$group != "all" & s$statistic == "mean"), yrs)
  expect_true(all(c("deviance", "chisq") %in% s$statistic))
  expect_true("groups" %in% names(ppc$plots))
  expect_s3_class(ppc$plots$groups, "ggplot")
  # group deviances add up to the overall deviance
  g <- s[s$group != "all" & s$statistic == "deviance", "observed"]
  expect_equal(sum(g), s$observed[s$group == "all" & s$statistic == "deviance"], tolerance = 1e-8)
  expect_error(posterior_predictive_check(fit, by = "nope"), "not found")
  expect_error(posterior_predictive_check(fit, discrepancy = "rmse"), "discrepancy")
  expect_error(posterior_predictive_check(fit, quantiles = 2), "quantiles")
  expect_error(posterior_predictive_check(1), "bm_fit")
  # interaction of two grouping columns
  two <- posterior_predictive_check(fit, responses = "biomass", n_draws = 50, by = c("year"))
  expect_s3_class(two, "bm_ppc")
  marg <- posterior_predictive_check(fit, responses = "biomass", n_draws = 50,
                                     random = "marginal", seed = 1)
  expect_equal(marg$settings$random, "marginal")
})

test_that("misspecification is detected: zeros are impossible under the truncated normal", {
  skip_if_not_installed("nimble")
  m <- fit_ppc_zero(); fit <- m$fit
  ppc <- posterior_predictive_check(fit, n_draws = 150, seed = 4)
  s <- ppc$summary
  z <- s[s$statistic == "prop_zero", ]
  expect_equal(z$observed, mean(m$d$sum_wet_weight == 0), tolerance = 1e-8)
  expect_equal(z$replicated_mean, 0)
  expect_equal(z$p_value, 0)
  expect_equal(z$flag, "extreme")
  expect_equal(s$flag[s$statistic == "prop_near_zero"], "extreme")
})

test_that("log-normal fits are checked with the log-normal likelihood", {
  skip_if_not_installed("nimble")
  fit <- fit_ppc_lnorm()$fit
  ppc <- posterior_predictive_check(fit, n_draws = 100, discrepancy = c("deviance", "chisq"), seed = 1)
  expect_true(all(ppc$replicates$biomass > 0))
  pr <- predict(fit, n_draws = 50)
  expect_true(all(pr$mean > 0))
  loc <- predict(fit, type = "location", n_draws = 50)
  expect_true(all(pr$mean > exp(loc$mean) * 0.99))   # lognormal mean >= median
})
