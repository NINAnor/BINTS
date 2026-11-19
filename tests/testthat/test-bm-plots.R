render_ok <- function(p) {
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  })
  print(p)
  TRUE
}
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
build <- function(p) ggplot2::ggplot_build(p)

test_that("the family cdf agrees with numerical integration", {
  f <- bm_families()$truncnorm
  m <- -1; s <- 2
  for (y in c(0, 0.5, 2, 6)) {
    num <- integrate(function(t) exp(f$loglik(t, m, s)), 0, y)$value
    expect_equal(f$cdf(y, m, s), num, tolerance = 1e-6)
  }
  y <- seq(0, 20, length.out = 50)
  cd <- f$cdf(y, -30, 2)
  expect_true(all(diff(cd) >= -1e-12) && all(cd >= 0 & cd <= 1))
  expect_equal(f$cdf(0, 1, 2), 0, tolerance = 1e-12)
  expect_equal(bm_families()$lognormal$cdf(2, 0, 1), stats::plnorm(2, 0, 1))
})

test_that("residuals() has the documented definitions and calibration", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  raw <- residuals(fit, n_draws = 100)
  pr <- predict(fit, n_draws = 100)
  expect_equal(raw$residual, raw$observed - raw$fitted)
  expect_equal(raw$fitted, pr$mean)
  expect_equal(raw$.row[raw$response == "biomass"], fit$prep$kept_rows)
  for (ty in c("standardised", "posterior_predictive")) {
    r <- residuals(fit, type = ty, n_draws = 150)
    expect_true(all(is.finite(r$residual)))
    for (rs in c("biomass", "richness")) {
      v <- r$residual[r$response == rs]
      expect_lt(abs(mean(v)), 0.15)
      expect_gt(stats::sd(v), 0.75)
      expect_lt(stats::sd(v), 1.25)
    }
  }
  expect_error(residuals(fit, type = "deviance"), "arg")
  one <- residuals(fit, responses = "richness", n_draws = 50)
  expect_equal(unique(one$response), "richness")
})

test_that("plot_observed_predicted shows the 1:1 line, metrics and optional intervals", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  p <- plot_observed_predicted(fit, "biomass", n_draws = 100)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("GeomAbline", "GeomPoint", "GeomText") %in% geoms(p)))
  expect_false("GeomLinerange" %in% geoms(p))
  n <- nrow(fit$prep$data)
  expect_equal(nrow(p$data), n)
  pr <- predict(fit, responses = "biomass", n_draws = 100)
  expect_equal(p$data$predicted, pr$mean)
  expect_equal(p$data$observed, fit$prep$y$biomass)
  e <- p$data$observed - p$data$predicted
  lab <- build(p)$data[[which(geoms(p) == "GeomText")]]$label
  expect_equal(lab, sprintf("RMSE = %s\nMAE = %s", format(sqrt(mean(e^2)), digits = 3),
                            format(mean(abs(e)), digits = 3)))
  expect_true(render_ok(p))
  pi <- plot_observed_predicted(fit, "biomass", interval = TRUE, n_draws = 100, metrics = FALSE)
  expect_true("GeomLinerange" %in% geoms(pi))
  expect_false("GeomText" %in% geoms(pi))
  # prediction intervals are wider than credible intervals
  pp <- plot_observed_predicted(fit, "biomass", type = "predictive", interval = TRUE,
                                n_draws = 100, seed = 1)
  expect_gt(mean(pp$data$upper - pp$data$lower), 2 * mean(pi$data$upper - pi$data$lower))
  pm <- plot_observed_predicted(fit, "biomass", statistic = "median", n_draws = 100)
  expect_equal(pm$data$predicted, pr$median)
})

test_that("plot_observed_predicted facets, square panels and validation", {
  skip_if_not_installed("nimble")
  m <- fit_joint(); fit <- m$fit
  p <- plot_observed_predicted(fit, c("biomass", "richness"), n_draws = 80)
  b <- build(p)
  expect_equal(nrow(b$layout$layout), 2L)
  for (pp in b$layout$panel_params) expect_equal(pp$x.range, pp$y.range, tolerance = 1e-8)
  py <- plot_observed_predicted(fit, "biomass", facet_by = "year", n_draws = 80)
  expect_equal(nrow(build(py)$layout$layout), length(unique(m$d$year)))
  pyr <- plot_observed_predicted(fit, c("biomass", "richness"), facet_by = "year", n_draws = 80)
  expect_equal(nrow(build(pyr)$layout$layout), 2L * length(unique(m$d$year)))
  expect_true(render_ok(pyr))
  expect_error(plot_observed_predicted(fit, "biomass", facet_by = "locality", max_facets = 10),
               "max_facets")
  expect_error(plot_observed_predicted(fit, "biomass", facet_by = "nope"), "facet_by")
  expect_error(plot_observed_predicted(fit, "abundance"), "fitted responses")
  expect_error(plot_observed_predicted(1, "biomass"), "bm_fit")
  # precomputed predictions are used as supplied
  pr <- predict(fit, n_draws = 50)
  pr$mean[pr$response == "biomass"] <- 7
  pc <- plot_observed_predicted(fit, "biomass", predictions = pr)
  expect_true(all(pc$data$predicted == 7))
  expect_error(plot_observed_predicted(fit, "biomass", predictions = pr[, 1:5]), "predict")
  nd <- predict(fit, newdata = data.frame(year = 2020, locality = "loc_001",
                                          sampling_name = "a-01", x1 = 0, x2 = 0), n_draws = 20)
  expect_error(plot_observed_predicted(fit, "biomass", predictions = nd), "observations used")
})

test_that("plot_residuals handles continuous, categorical and mixed covariates", {
  skip_if_not_installed("nimble")
  m <- fit_joint(); fit <- m$fit
  pc <- plot_residuals(fit, c("longitude", "latitude", "x1"), "biomass", n_draws = 80)
  expect_s3_class(pc, "ggplot")
  expect_true(all(c("GeomHline", "GeomSmooth", "GeomPoint") %in% geoms(pc)))
  expect_equal(unique(build(pc)$data[[which(geoms(pc) == "GeomHline")]]$yintercept), 0)
  expect_equal(nrow(build(pc)$layout$layout), 3L)
  expect_true(render_ok(pc))
  nosm <- plot_residuals(fit, "x1", "biomass", smooth = FALSE, n_draws = 80)
  expect_false("GeomSmooth" %in% geoms(nosm))
  # categorical only: genuine discrete x scale, boxplots for few levels
  pk <- plot_residuals(fit, c("year", "region_name_absent")[1], "biomass", n_draws = 80)
  expect_true("GeomBoxplot" %in% geoms(pk))
  expect_false("GeomSmooth" %in% geoms(pk))
  expect_s3_class(build(pk)$layout$panel_scales_x[[1]], "ScaleDiscretePosition")
  pl <- plot_residuals(fit, "locality", "biomass", n_draws = 80, max_levels = 10)
  expect_false("GeomBoxplot" %in% geoms(pl))
  expect_true(render_ok(pl))
  # mixed: categorical at integer positions, levels in the facet label
  pm <- plot_residuals(fit, c("x1", "year"), c("biomass", "richness"), n_draws = 80)
  expect_equal(nrow(build(pm)$layout$layout), 4L)
  labs <- levels(pm$data$.covariate)
  expect_true(any(grepl("year \\(5 levels\\)", labs)))
  expect_true(render_ok(pm))
  # residual flavours
  for (ty in c("standardised", "posterior_predictive")) {
    pt <- plot_residuals(fit, "x1", "biomass", type = ty, n_draws = 80)
    expect_match(pt$labels$y, if (ty == "standardised") "Standardised" else "normal scores")
  }
  expect_error(plot_residuals(fit, "nope", "biomass"), "not found")
  expect_error(plot_residuals(fit, character(), "biomass"), "covariates")
  expect_error(plot_residuals(fit, "x1", "abundance"), "fitted responses")
})

test_that("plot_response_covariates compares observed and predicted on common axes", {
  skip_if_not_installed("nimble")
  m <- fit_joint(); fit <- m$fit
  p <- plot_response_covariates(fit, c("x1", "x2"), c("biomass", "richness"), n_draws = 80)
  expect_s3_class(p, "ggplot")
  b <- build(p)
  lay <- b$layout$layout
  expect_equal(nrow(lay), 2L * 2L * 2L)
  expect_setequal(levels(p$data$.source), c("Observed", "Predicted"))
  # y limits of observed and predicted panels are equal within a response
  yr <- vapply(b$layout$panel_params, function(pp) paste(round(pp$y.range, 6), collapse = "~"), "")
  for (r in c("biomass", "richness")) {
    ii <- which(lay$response == r)
    expect_equal(length(unique(yr[ii])), 1L)
  }
  expect_true("GeomSmooth" %in% geoms(p))
  expect_false("GeomSmooth" %in% geoms(plot_response_covariates(fit, "x1", "biomass", trend = FALSE,
                                                                 n_draws = 50)))
  # predicted values are the fitted expected responses, observed the data
  one <- plot_response_covariates(fit, "x1", "biomass", n_draws = 80)
  d <- one$data
  expect_equal(d$.y[d$.source == "Observed"], fit$prep$y$biomass)
  expect_equal(d$.y[d$.source == "Predicted"], predict(fit, responses = "biomass", n_draws = 80)$mean)
  expect_true(render_ok(p))
  pm <- plot_response_covariates(fit, c("x1", "year"), "biomass", n_draws = 50)
  expect_true(render_ok(pm))
})

test_that("covariate effect curves are model-based and distinct from fitted values", {
  skip_if_not_installed("nimble")
  fit <- fit_joint()$fit
  for (ef in c("conditional", "marginal")) {
    p <- plot_response_covariates(fit, "x1", "biomass", effect = ef, n_draws = 120,
                                  n_rows = 80, n_grid = 12, seed = 1)
    expect_true(all(c("GeomRibbon", "GeomLine") %in% geoms(p)))
    ln <- build(p)$data[[which(geoms(p) == "GeomLine")]]
    expect_equal(nrow(ln), 12L)
    expect_gt(ln$y[12], ln$y[1] + 3)    # biomass increases with x1 (true effect 5 per sd)
    rb <- build(p)$data[[which(geoms(p) == "GeomRibbon")]]
    expect_true(all(rb$ymin <= ln$y & ln$y <= rb$ymax))
    expect_match(p$labels$caption, ef)
    expect_true(render_ok(p))
  }
  # richness has no x1 effect in truth; its x2 effect is positive
  p2 <- plot_response_covariates(fit, "x2", "richness", effect = "marginal", n_draws = 120,
                                 n_rows = 80, n_grid = 10, seed = 1)
  ln <- build(p2)$data[[which(geoms(p2) == "GeomLine")]]
  expect_gt(ln$y[10], ln$y[1])
  # only model covariates have effect curves
  expect_error(plot_response_covariates(fit, "latitude", "biomass", effect = "marginal"),
               "Effect curves")
  expect_error(plot_response_covariates(fit, "year", "biomass", effect = "conditional"),
               "Effect curves")
})

test_that("plot_spatial_predictions maps precomputed predictions without refitting", {
  skip_if_not_installed("nimble")
  m <- fit_spatial(); fit <- m$fit
  bb <- sf::st_bbox(c(xmin = 5, xmax = 15, ymin = 58, ymax = 68), crs = 4326)
  cells <- sf::st_make_grid(sf::st_as_sfc(bb), n = c(6, 6))
  grid <- do.call(rbind, lapply(2021:2022, function(y) sf::st_sf(year = y, x1 = 0, geometry = cells)))
  pred <- predict(fit, grid, random = "zero", n_draws = 60)
  expect_s3_class(pred, "sf")
  expect_no_message(p <- plot_spatial_predictions(fit, pred, "biomass"))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 36L * 2L)
  expect_equal(p$data$.value, p$data$median)
  expect_equal(nrow(build(p)$layout$layout), 2L)      # one panel per year
  expect_true(render_ok(p))
  # statistics and uncertainty displays
  expect_equal(plot_spatial_predictions(fit, pred, "biomass", statistic = "mean")$data$.value,
               pred$mean)
  expect_equal(plot_spatial_predictions(fit, pred, "biomass", uncertainty = TRUE)$data$.value,
               pred$upper - pred$lower)
  expect_equal(plot_spatial_predictions(fit, pred, "biomass", uncertainty = "sd")$data$.value,
               pred$sd)
  pa <- plot_spatial_predictions(fit, pred, "biomass", uncertainty = "alpha")
  expect_true(render_ok(pa))
  expect_error(plot_spatial_predictions(fit, pred, "biomass", uncertainty = "var"), "uncertainty")
  # year subset, boundary and display CRS
  p1 <- plot_spatial_predictions(fit, pred, "biomass", year = 2022)
  expect_equal(nrow(p1$data), 36L)
  expect_equal(nrow(build(p1)$layout$layout), 1L)
  expect_error(plot_spatial_predictions(fit, pred, "biomass", year = 1999), "No prediction rows")
  bd <- sf::st_transform(sf::st_as_sf(sf::st_as_sfc(bb)), 25833)
  pb <- plot_spatial_predictions(fit, pred, "biomass", boundary = bd, crs = 25833)
  expect_equal(length(pb$layers), length(p$layers) + 1L)
  expect_true(render_ok(pb))
  expect_error(plot_spatial_predictions(fit, pred, "biomass", boundary = 1), "boundary")
  expect_error(plot_spatial_predictions(fit, pred, "richness"), "fitted responses")
  expect_error(plot_spatial_predictions(fit, pred, "biomass", point_size = -1), "point_size")
  expect_error(plot_spatial_predictions(fit, 1:3, "biomass"), "sf object or a data frame")
  # an ordinary ggplot2 object: users can keep customising it
  expect_s3_class(p + ggplot2::ggtitle("x") + ggplot2::theme_bw(), "ggplot")
})

test_that("plot_spatial_predictions computes predictions only when asked to", {
  skip_if_not_installed("nimble")
  fit <- fit_spatial()$fit
  bb <- sf::st_bbox(c(xmin = 6, xmax = 14, ymin = 59, ymax = 67), crs = 4326)
  cells <- sf::st_make_grid(sf::st_as_sfc(bb), n = c(4, 4))
  grid <- sf::st_sf(year = 2022, x1 = 0, geometry = cells)
  expect_message(p <- plot_spatial_predictions(fit, grid, "biomass",
                                               predict_args = list(n_draws = 40)),
                 "Computing predictions")
  expect_equal(nrow(p$data), 16L)
  expect_true(all(p$data$mean > 0))
  # points as a data frame of coordinates (the model's coordinate columns)
  df <- data.frame(longitude = c(7, 9, 11, 13), latitude = c(60, 62, 64, 66), year = 2021, x1 = 0)
  expect_message(pp <- plot_spatial_predictions(fit, df, "biomass",
                                                predict_args = list(n_draws = 40)), "Computing")
  expect_true(all(sf::st_geometry_type(pp$data) == "POINT"))
  expect_true(render_ok(pp))
  expect_error(plot_spatial_predictions(fit, data.frame(a = 1), "biomass"), "x1|coordinate")
})

test_that("several responses can share a map with relative scaling", {
  skip_if_not_installed("nimble")
  fit <- fit_spatial()$fit
  cells <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(c(xmin = 5, xmax = 15, ymin = 58, ymax = 68),
                                                      crs = 4326)), n = c(5, 5))
  base <- sf::st_sf(year = rep(c(2021, 2022), each = 25), geometry = rep(cells, 2))
  ctr <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(base)))
  mk <- function(r, s) {
    d <- base
    d$response <- r
    d$median <- s * (20 + ctr[, 2] - 58)
    d$mean <- d$median
    d$sd <- s
    d$lower <- d$median - 3 * s
    d$upper <- d$median + 3 * s
    d
  }
  # a response that was not fitted cannot be requested
  pd <- rbind(mk("biomass", 1), mk("biomass", 1))
  pd2 <- pd
  pd2$response <- "richness"
  # (use biomass for both rows of data, relabelled, to test the layout code)
  both <- rbind(mk("biomass", 1), transform(mk("biomass", 25), response = "richness"))
  fit2 <- fit
  fit2$spec$responses <- c("biomass", "richness")
  p <- plot_spatial_predictions(fit2, both, c("biomass", "richness"), scale = "relative")
  expect_true(all(p$data$.value >= 0 & p$data$.value <= 1))
  for (r in c("biomass", "richness")) {
    v <- p$data$.value[p$data$response == r]
    expect_equal(range(v), c(0, 1))
  }
  expect_equal(nrow(build(p)$layout$layout), 4L)
  pr <- plot_spatial_predictions(fit2, both, c("biomass", "richness"), facet_by = "response")
  expect_equal(nrow(build(pr)$layout$layout), 4L)
  expect_true(render_ok(p))
})
