test_that("count families have consistent pmf, moments, cdf and simulation", {
  for (fam in c("poisson", "negbinomial")) {
    f <- bm_families()[[fam]]
    for (m in log(c(0.5, 6, 300))) {
      y <- matrix(0:4000, ncol = 1)
      mu <- matrix(m, nrow(y), 1)
      p <- exp(f$loglik(y, mu, 3))
      expect_equal(sum(p), 1, tolerance = 1e-8)
      expect_equal(sum(y * p), f$expected(matrix(m), 3)[1], tolerance = 1e-6)
      expect_equal(sum((y - exp(m))^2 * p), f$variance(matrix(m), 3)[1], tolerance = 1e-5)
      # cdf - left limit is the pmf
      expect_equal(as.vector(f$cdf(y, mu, 3) - f$cdf_left(y, mu, 3)), as.vector(p), tolerance = 1e-10)
    }
    set.seed(3)
    s <- f$simulate(matrix(log(40), 1e5, 1), 3)
    expect_true(is.matrix(s) && all(s == round(s)) && all(s >= 0))
    expect_equal(mean(s), 40, tolerance = 0.02)
    expect_equal(stats::var(as.vector(s)), f$variance(matrix(log(40)), 3)[1], tolerance = 0.05)
  }
  nb <- bm_families()$negbinomial
  po <- bm_families()$poisson
  # large size: the negative binomial approaches the Poisson
  expect_equal(nb$loglik(5, log(8), 1e8), po$loglik(5, log(8), NA), tolerance = 1e-5)
  # overdispersion: smaller size, larger variance
  expect_gt(nb$variance(log(50), 2), nb$variance(log(50), 20))
  expect_equal(nb$variance(log(50), 1e9), po$variance(log(50), NA), tolerance = 1e-4)
})

test_that("each response accepts only its families", {
  both <- c("biomass", "richness")
  ok <- function(f) bm_resolve_family(f, both)
  expect_equal(ok(c(biomass = "lognormal", richness = "poisson"))[["richness"]], "poisson")
  expect_equal(ok(c(biomass = "truncnorm", richness = "negbinomial"))[["richness"]], "negbinomial")
  expect_equal(ok(c(biomass = "truncnorm", richness = "truncnorm"))[["richness"]], "truncnorm")
  expect_error(ok(c(biomass = "poisson", richness = "poisson")), "not available for response .biomass.")
  expect_error(ok(c(biomass = "negbinomial", richness = "negbinomial")), "use one of")
  expect_error(ok(c(biomass = "truncnorm", richness = "lognormal")), "not available for response .richness.")
  expect_error(ok("gamma"), "Unknown family")
  expect_error(fit_model(bm_test_data(), responses = "biomass", family = "poisson"), "not available")
})

test_that("generated code uses the right distribution and parameters per family", {
  d <- bm_nb_data()
  mk <- function(fam, fixed = list(richness = "x2")) {
    spec <- list(responses = "richness", family = c(richness = fam), fixed_effects = fixed,
                 random_effects = "year", temporal_effect = "none", spatial_effect = "none")
    prep <- bm_prepare_data(d, "richness", fixed, "year", spec$family, bm_default_columns())
    pri <- bm_default_priors(prep, spec$family)
    list(b = bm_build_model(prep, spec, pri), prep = prep, spec = spec, pri = pri)
  }
  txt <- function(m) gsub("\\s+", " ", paste(deparse(m$b$code), collapse = " "))
  nb <- mk("negbinomial")
  expect_match(txt(nb), "prob_richness\\[i\\] <- phi_richness/\\(phi_richness \\+ exp\\(mu_richness\\[i\\]\\)\\)")
  expect_match(txt(nb), "richness\\[i\\] ~ dnbinom\\(prob = prob_richness\\[i\\], size = phi_richness\\)")
  expect_match(txt(nb), "phi_richness ~ dgamma\\(prior_disp_shape_richness, prior_disp_rate_richness\\)")
  expect_true("phi_richness" %in% nb$b$monitors)
  expect_false(any(grepl("sigma", nb$b$monitors)))
  expect_equal(nb$b$constants$prior_disp_shape_richness, 1)
  expect_equal(nb$b$constants$prior_disp_rate_richness, 0.1)
  # the log mean count has a sensible prior location
  expect_equal(nb$pri$richness$beta0_mean, mean(log1p(d$no_species)))
  expect_true(all(c("phi_richness") %in% names(bm_make_inits(nb$prep, nb$spec, nb$pri, 1))))
  # the constant p_richness (number of covariates) must not collide with a node
  expect_false("p_richness" %in% all.vars(nb$b$code) && "p_richness" %in% names(nb$b$data))
  po <- mk("poisson")
  expect_match(txt(po), "lambda_richness\\[i\\] <- exp\\(mu_richness\\[i\\]\\)")
  expect_match(txt(po), "richness\\[i\\] ~ dpois\\(lambda_richness\\[i\\]\\)")
  expect_false(grepl("phi_|sigma_", txt(po)))
  expect_false(any(grepl("phi|sigma", po$b$monitors)))
  expect_false(any(grepl("phi|sigma", names(bm_make_inits(po$prep, po$spec, po$pri, 1)))))
})

test_that("count families require non-negative integers; priors can be overridden", {
  d <- bm_nb_data()
  d$no_species[3] <- 10.5
  expect_error(fit_model(d, responses = "richness", family = "poisson"), "non-negative integer counts")
  expect_error(fit_model(d, responses = "richness", family = "negbinomial"), "non-negative integer counts")
  d$no_species[3] <- -1
  expect_error(fit_model(d, responses = "richness", family = "negbinomial"), "integer counts")
  d <- bm_nb_data()
  pr <- bm_resolve_priors(list(richness = list(disp_shape = 2, disp_rate = 0.5)),
                          bm_default_priors(bm_prepare_data(d, "richness", list(richness = character()),
                            "year", c(richness = "negbinomial"), bm_default_columns()), c(richness = "negbinomial")))
  expect_equal(unname(unlist(pr$richness[c("disp_shape", "disp_rate")])), c(2, 0.5))
  expect_error(bm_resolve_priors(list(richness = list(disp_shape = -1)), pr), "scales positive")
})

test_that("negative binomial and Poisson models fit, predict and residualise counts", {
  skip_if_not_installed("nimble")
  m <- fit_nb(); fit <- m$fit
  expect_true(all(c("phi_richness", "beta0_richness", "sd_locality_richness") %in% fit$diagnostics$parameter))
  expect_output(print(fit), "negbinomial")
  loc <- predict(fit, type = "location", summarise = FALSE, n_draws = 60)
  ex <- predict(fit, type = "expected", summarise = FALSE, n_draws = 60)
  expect_equal(ex$draws$richness, exp(loc$draws$richness))
  pr <- predict(fit, type = "predictive", summarise = FALSE, n_draws = 60, seed = 1)
  expect_true(all(pr$draws$richness == round(pr$draws$richness)) && all(pr$draws$richness >= 0))
  expect_equal(mean(pr$draws$richness), mean(ex$draws$richness), tolerance = 0.03)
  # predictive intervals are wider than intervals for the expectation
  pe <- predict(fit, n_draws = 100); pp <- predict(fit, type = "predictive", n_draws = 100, seed = 2)
  expect_gt(mean(pp$upper - pp$lower), 2 * mean(pe$upper - pe$lower))
  nd <- data.frame(year = 2999, locality = "new", sampling_name = "a-01", x2 = c(-1, 1))
  pn <- predict(fit, nd, random = "marginal", n_draws = 100, seed = 1)
  expect_gt(pn$mean[2], pn$mean[1])      # positive true covariate effect
  r <- residuals(fit, type = "posterior_predictive", n_draws = 100)
  expect_lt(abs(mean(r$residual)), 0.2)
  expect_gt(stats::sd(r$residual), 0.75); expect_lt(stats::sd(r$residual), 1.25)
  expect_s3_class(plot_residuals(fit, "x2", "richness", type = "posterior_predictive", n_draws = 60), "ggplot")
  expect_s3_class(plot_observed_predicted(fit, "richness", n_draws = 60), "ggplot")
  # Poisson has no second parameter
  pf <- fit_pois()$fit
  expect_false(any(grepl("^(phi|sigma)", pf$diagnostics$parameter)))
  pp <- predict(pf, n_draws = 50)
  expect_true(all(pp$mean > 0))
})

test_that("the posterior predictive check detects overdispersion that a Poisson model ignores", {
  skip_if_not_installed("nimble")
  ppo <- posterior_predictive_check(fit_pois()$fit, n_draws = 200, seed = 1)$summary
  ppn <- posterior_predictive_check(fit_nb()$fit, n_draws = 200, seed = 1)$summary
  get <- function(s, st) s[s$statistic == st & s$group == "all", ]
  expect_equal(get(ppo, "deviance")$flag, "extreme")
  expect_equal(get(ppo, "sd")$flag, "extreme")
  expect_gt(get(ppo, "sd")$observed, get(ppo, "sd")$replicated_upper)
  expect_equal(get(ppn, "deviance")$flag, "")
  expect_equal(get(ppn, "sd")$flag, "")
  # replicates are counts
  rep <- posterior_predictive_check(fit_nb()$fit, n_draws = 50, n_overlay = 20, seed = 1)$replicates$richness
  expect_true(all(rep == round(rep)))
})

test_that("a joint model combines a log-normal biomass with a negative binomial richness", {
  skip_if_not_installed("nimble")
  fit <- fit_counts()$fit
  expect_equal(unname(fit$spec$family), c("lognormal", "negbinomial"))
  expect_true(all(c("sigma_biomass", "phi_richness") %in% fit$diagnostics$parameter))
  p <- predict(fit, n_draws = 60)
  expect_true(all(p$mean > 0))
  ppc <- posterior_predictive_check(fit, n_draws = 60, discrepancy = c("deviance", "chisq"), seed = 1)
  expect_setequal(unique(ppc$summary$response), c("biomass", "richness"))
})

test_that("negative binomial parameters are recovered (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(
    n_year = 5, n_locality = 60, n_rep = 4, seed = 41,
    family = c(biomass = "truncnorm", richness = "negbinomial"),
    beta0 = c(biomass = 20, richness = log(500)),
    beta = list(biomass = c(x1 = 5, x2 = 0), richness = c(x1 = 0.2, x2 = -0.3)),
    sd_year = c(biomass = 3, richness = 0.1), sd_locality = c(biomass = 6, richness = 0.3),
    sd_replication = c(biomass = 4, richness = 0.1), nb_size = 6)
  fit <- bm_test_fit(d, responses = "richness", family = "negbinomial", fixed_effects = c("x1", "x2"),
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  est <- fit$diagnostics
  cov <- function(p, v) { r <- est[est$parameter == p, ]; r$q2.5 <= v && v <= r$q97.5 }
  expect_true(cov("phi_richness", 6))
  expect_true(cov("beta_richness[1]", 0.2))
  expect_true(cov("beta_richness[2]", -0.3))
  expect_true(cov("sd_locality_richness", 0.3))
})

