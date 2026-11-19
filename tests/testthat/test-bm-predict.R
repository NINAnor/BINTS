# Fitted models shared by the prediction tests (fitted once)
test_that("predict() returns tidy posterior summaries for the observations", {
  skip_if_not_installed("nimble")
  m <- fit_joint(); fit <- m$fit
  p <- predict(fit, n_draws = 200)
  n <- nrow(fit$prep$data)
  expect_equal(nrow(p), 2L * n)
  expect_true(all(c(".row", "response", "mean", "median", "sd", "lower", "upper") %in% names(p)))
  expect_setequal(unique(p$response), c("biomass", "richness"))
  expect_equal(p$.row[p$response == "biomass"], fit$prep$kept_rows)
  expect_true(all(p$lower <= p$median & p$median <= p$upper))
  expect_true(all(p$mean > 0))
  one <- predict(fit, responses = "richness", n_draws = 50)
  expect_equal(unique(one$response), "richness")
  expect_error(predict(fit, responses = "foo"), "responses")
  expect_error(predict(fit, level = 1.2), "level")
  wide <- predict(fit, n_draws = 200, level = 0.99)
  expect_true(all(wide$upper - wide$lower >= p$upper - p$lower - 1e-8))
})

test_that("expected response is the truncated mean of the location draws", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  loc <- predict(fit, type = "location", summarise = FALSE, n_draws = 60)
  exp_ <- predict(fit, type = "expected", summarise = FALSE, n_draws = 60)
  expect_s3_class(loc, "bm_draws")
  mat <- as.matrix(fit$samples)[loc$draw_index, ]
  for (r in c("biomass", "richness")) {
    manual <- bm_families()$truncnorm$expected(loc$draws[[r]], mat[, paste0("sigma_", r)])
    expect_equal(exp_$draws[[r]], manual)
    expect_true(all(exp_$draws[[r]] >= loc$draws[[r]]))
    expect_true(any(exp_$draws[[r]] > loc$draws[[r]]))
    expect_true(all(exp_$draws[[r]] > 0))
  }
})

test_that("predictive draws are non-negative and centred on the expectation", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  pr <- predict(fit, type = "predictive", summarise = FALSE, n_draws = 100, seed = 3)
  ex <- predict(fit, type = "expected", summarise = FALSE, n_draws = 100)
  for (r in c("biomass", "richness")) {
    expect_true(all(pr$draws[[r]] >= 0))
    expect_equal(mean(pr$draws[[r]]), mean(ex$draws[[r]]), tolerance = 0.02)
    # observation noise makes predictive draws more spread out than expectations
    expect_gt(stats::sd(pr$draws[[r]][, 1]), stats::sd(ex$draws[[r]][, 1]))
  }
  pr2 <- predict(fit, type = "predictive", summarise = FALSE, n_draws = 100, seed = 3)
  expect_equal(pr$draws, pr2$draws)
})

test_that("newdata with known levels reproduces the fitted predictions", {
  skip_if_not_installed("nimble")
  m <- fit_joint(); fit <- m$fit
  base <- predict(fit, n_draws = 100)
  nd <- predict(fit, newdata = fit$prep$data[1:12, ], n_draws = 100)
  b <- base[base$.row %in% 1:12, c(".row", "response", "mean", "lower", "upper")]
  n2 <- nd[, c(".row", "response", "mean", "lower", "upper")]
  expect_equal(n2$mean, b$mean)
  expect_equal(n2$upper, b$upper)
})

test_that("missing and invalid newdata are handled", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  nd <- fit$prep$data[1:5, ]
  expect_error(predict(fit, nd[, setdiff(names(nd), "x1")]), "x1")
  nd$x1[2] <- NA
  expect_message(p <- predict(fit, nd, responses = "biomass", n_draws = 50), "missing")
  expect_true(is.na(p$mean[2]) && !any(is.na(p$mean[-2])))
  expect_error(predict(fit, 1:3), "newdata")
  expect_error(predict(fit, nd[0, ]), "no rows")
})

test_that("random-effect modes and new levels behave as documented", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  nd <- data.frame(year = c(2020, 2999), locality = c("loc_001", "brand_new"),
                   sampling_name = "x_sampling-01", x1 = 0, x2 = 0)
  get <- function(...) predict(fit, nd, responses = "biomass", n_draws = 150, ...)
  expect_error(get(new_levels = "error"), "not in the fitting data")
  z1 <- get(new_levels = "zero"); z2 <- get(new_levels = "zero")
  expect_equal(z1$mean, z2$mean)
  s1 <- get(seed = 4); s2 <- get(seed = 4); s3 <- get(seed = 5)
  expect_equal(s1$mean, s2$mean)
  expect_false(isTRUE(all.equal(s1$mean, s3$mean)))
  # simulated new levels carry extra uncertainty relative to a known level
  expect_gt(s1$sd[2], s1$sd[1])
  expect_gt(s1$sd[2], z1$sd[2])
  # zero-effect mode ignores all random effects
  zr <- get(random = "zero")
  expect_equal(zr$mean[1], zr$mean[2])
  # named modes: only year is marginalised, locality stays conditional
  mix <- get(random = c(year = "zero"), seed = 1)
  expect_equal(mix$mean[1], get(random = c(year = "zero"), seed = 2)$mean[1])
  expect_error(get(random = c(site = "zero")), "random effects")
  expect_error(get(random = "average"), "random")
  expect_error(get(random = c("zero", "marginal")), "single mode")
  # marginal: reproducible, and wider than the all-zero typical group
  m1 <- get(random = "marginal", seed = 8); m2 <- get(random = "marginal", seed = 8)
  expect_equal(m1$mean, m2$mean)
  expect_gt(m1$upper[1] - m1$lower[1], zr$upper[1] - zr$lower[1])
})

test_that("new random-effect draws follow shared and correlated structure", {
  mk <- function(sharing) {
    D <- 1
    mat <- matrix(c(2, 5, 3, 10, 0.7, 4), nrow = D,
                  dimnames = list(NULL, c("sd_year_biomass", "sd_year_richness",
                                          "lambda_re_year_richness", "sd_locality_biomass",
                                          "rho_locality_richness", "sd_locality_richness")))
    mat <- mat[rep(1, 1), , drop = FALSE]
    list(fit = list(spec = list(responses = c("biomass", "richness"), sharing = sharing)),
         mat = mat)
  }
  sh <- list(active = TRUE, loadings = TRUE, random = "year", correlated = "locality")
  z <- mk(sh)
  set.seed(1)
  out <- bm_sim_new_effects(z$fit, "year", 1000, z$mat)
  expect_equal(out$richness, 3 * out$biomass)
  expect_equal(stats::sd(out$biomass), 2, tolerance = 0.1)
  out2 <- bm_sim_new_effects(z$fit, "locality", 40000, z$mat)
  expect_equal(stats::cor(out2$biomass[1, ], out2$richness[1, ]), 0.7, tolerance = 0.02)
  expect_equal(stats::sd(out2$richness), 4, tolerance = 0.1)
  sh2 <- list(active = TRUE, loadings = FALSE, random = "year", correlated = character())
  out3 <- bm_sim_new_effects(mk(sh2)$fit, "year", 10, z$mat)
  expect_equal(out3$richness, out3$biomass)
})

test_that("spatial prediction works at known sites, new points and sf geometries", {
  skip_if_not_installed("nimble")
  m <- fit_spatial(); fit <- m$fit; d <- m$d
  site <- unique(d[, c("locality", "longitude", "latitude")])[1:3, ]
  nd_known <- data.frame(site, year = 2022, x1 = 0)
  nd_xy <- data.frame(longitude = site$longitude, latitude = site$latitude, year = 2022, x1 = 0)
  kn <- predict(fit, nd_known, random = "zero", n_draws = 100)
  xy <- predict(fit, nd_xy, random = "zero", n_draws = 100)
  # the spline is deterministic: known site and the same coordinates agree
  expect_equal(kn$mean, xy$mean)
  # new points, including one outside the sampled area: finite, positive
  nd_new <- data.frame(longitude = c(7, 12, 30), latitude = c(60, 66, 80), year = 2023, x1 = 0)
  pn <- predict(fit, nd_new, random = "zero", n_draws = 100)
  expect_true(all(is.finite(pn$mean)) && all(pn$mean > 0))
  # sf input (lon/lat and another CRS) gives the same predictions
  pts <- sf::st_as_sf(data.frame(nd_new[, c("year", "x1")], g = 1:3,
                                 longitude = nd_new$longitude, latitude = nd_new$latitude),
                      coords = c("longitude", "latitude"), crs = 4326)
  ps <- predict(fit, pts, random = "zero", n_draws = 100)
  expect_s3_class(ps, "sf")
  expect_equal(ps$mean, pn$mean, tolerance = 1e-6)
  ps2 <- predict(fit, sf::st_transform(pts, 25833), random = "zero", n_draws = 100)
  expect_equal(ps2$mean, pn$mean, tolerance = 1e-4)
  # coordinates needed for unknown locations
  expect_error(predict(fit, data.frame(year = 2022, x1 = 0)), "coordinate")
  # spatial = "none" drops the field
  p0 <- predict(fit, nd_known, random = "zero", spatial = "none", n_draws = 100)
  expect_false(isTRUE(all.equal(p0$mean, kn$mean)))
})

test_that("a linear trend extrapolates to new years", {
  skip_if_not_installed("nimble")
  m <- fit_spatial(); fit <- m$fit; d <- m$d
  site <- unique(d[, c("locality", "longitude", "latitude")])[1, ]
  nd <- data.frame(site, year = c(2025, 2035), x1 = 0, row.names = NULL)
  lc <- predict(fit, nd, type = "location", random = "zero", summarise = FALSE, n_draws = 80)
  beta <- as.matrix(fit$samples)[lc$draw_index, "beta_trend_biomass"]
  expect_equal(lc$draws$biomass[, 2] - lc$draws$biomass[, 1], 10 * beta, tolerance = 1e-8)
  expect_error(predict(fit, data.frame(site, x1 = 0)), "year")
})

test_that("cyclic effect is periodic and range-checked", {
  skip_if_not_installed("nimble")
  fit <- fit_cyclic()$fit
  nd <- data.frame(year = 2021, locality = "loc_001", sampling_name = "x_sampling-01",
                   start_date_julian = c(0, 365, 150))
  lc <- predict(fit, nd, type = "location", random = "zero", summarise = FALSE, n_draws = 50)
  expect_equal(lc$draws$biomass[, 1], lc$draws$biomass[, 2], tolerance = 1e-8)
  expect_false(isTRUE(all.equal(lc$draws$biomass[, 1], lc$draws$biomass[, 3])))
  nd$start_date_julian[3] <- 400
  expect_error(predict(fit, nd), "within")
})

test_that("NNGP predictions at unsampled locations include spatial uncertainty", {
  skip_if_not_installed("nimble")
  m <- fit_nngp(); fit <- m$fit
  nd <- data.frame(longitude = c(8, 8), latitude = c(62, 62), year = 2021)
  e1 <- predict(fit, nd, random = "zero", n_draws = 150, seed = 1)
  expect_equal(e1$mean[1], e1$mean[2])     # same coordinates share one draw
  site <- unique(m$d[, c("locality", "longitude", "latitude")])[1, ]
  known <- predict(fit, data.frame(site, year = 2021), random = "zero", n_draws = 150)
  newpt <- predict(fit, data.frame(site[, 2:3] + c(0.5, 0.5), year = 2021), random = "zero",
                   n_draws = 150, seed = 2)
  expect_gt(newpt$sd, known$sd)
  expect_true(all(is.finite(e1$mean)) && all(e1$mean > 0))
})
