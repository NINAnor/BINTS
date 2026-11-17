sp_data <- function(type, n_locality = 30, seed = 55, ...) {
  simulate_bm_data(
    n_year = 3, n_locality = n_locality, n_rep = 3, seed = seed, spatial = type,
    spatial_sd = c(biomass = 8, richness = 100), spatial_range_km = 250,
    sd_locality = c(biomass = 0.01, richness = 1),
    sd_replication = c(biomass = 0.01, richness = 1), ...
  )
}

test_that("spatial terms enter the model code only when requested", {
  d <- sp_data("autocorrelation")
  fam <- c(biomass = "truncnorm")
  mk <- function(se) {
    spec <- list(responses = "biomass", family = fam,
                 fixed_effects = list(biomass = character()),
                 random_effects = "year", temporal_effect = "none",
                 spatial_effect = se)
    prep <- bm_prepare_data(d, "biomass", spec$fixed_effects, "year", fam,
                            bm_default_columns(), spatial_effect = se)
    bm_build_model(prep, spec, bm_default_priors(prep, fam))
  }
  txt <- function(b) gsub("\\s+", " ", paste(deparse(b$code), collapse = " "))
  expect_false(grepl("w_sp|f_sp", txt(mk("none"))))
  ac <- mk("autocorrelation")
  expect_match(txt(ac), "phi_idx_sp_biomass ~ dcat")
  expect_match(txt(ac), "w_sp_biomass\\[site_id\\[i\\]\\]")
  expect_true(all(c("Bnn", "sqrtF") %in% names(ac$data)))
  co <- mk("coordinates")
  expect_match(txt(co), "Bsp\\[1:L_site, 1:K_sp\\] %\\*% v_sp_biomass")
  expect_equal(nrow(co$constants$Bsp), co$constants$L_site)
})

test_that("NNGP coefficients equal the exact Gaussian-process conditionals", {
  d <- sp_data("autocorrelation")
  prep <- bm_prepare_data(d, "biomass", list(biomass = character()), "year",
                          c(biomass = "truncnorm"), bm_default_columns(),
                          spatial_effect = "autocorrelation")
  ng <- prep$spatial$nngp
  xy <- ng$xy
  expect_true(all(ng$nbr[2:6, 1] < 2:6))
  # sites 2..6 have fewer than m earlier sites: neighbour set is all earlier sites
  g <- 8
  ell <- ng$range_grid[g]
  D <- bm_distance_km(xy, xy, TRUE)
  C <- exp(-D / ell)
  i <- 6
  k <- i - 1
  b_exact <- solve(C[1:k, 1:k] + diag(1e-6, k), C[1:k, i])
  b <- numeric(k)
  b[ng$nbr[i, 1:k]] <- ng$Bnn[i, 1:k, g]
  expect_equal(b, b_exact, tolerance = 1e-6)
  f_exact <- 1 - sum(C[1:k, i] * b_exact)
  expect_equal(ng$sqrtF[i, g]^2, f_exact, tolerance = 1e-6)
  expect_true(all(ng$sqrtF > 0 & ng$sqrtF <= 1))
  expect_true(all(diff(ng$range_grid) > 0))
})

test_that("coordinate spline basis is centred and predicts at new points", {
  d <- sp_data("coordinates", n_locality = 40)
  prep <- bm_prepare_data(d, "biomass", list(biomass = character()), "year",
                          c(biomass = "truncnorm"), bm_default_columns(),
                          spatial_effect = "coordinates")
  tp <- prep$spatial$tp
  expect_equal(unname(colSums(tp$basis)), rep(0, tp$K), tolerance = 1e-6)
  expect_gte(tp$K - tp$Kp, 2L)
  P <- bm_plane(prep$spatial$xy, tp$info)
  B <- mgcv::PredictMat(tp$smooth, data.frame(cx = P[, 1], cy = P[, 2])) %*% tp$transform
  expect_equal(unname(B), unname(tp$basis), tolerance = 1e-8)
})

test_that("spatial inputs are validated", {
  d <- sp_data("autocorrelation")
  expect_error(fit_model(d, responses = "biomass", spatial_effect = "autocorrelation",
                         coordinates = "longitude"), "two column names")
  expect_error(fit_model(d, responses = "biomass", spatial_effect = "autocorrelation",
                         coordinates = c("lon", "lat")), "not found")
  d2 <- d
  d2$longitude[1] <- d2$longitude[1] + 1
  expect_error(fit_model(d2, responses = "biomass", spatial_effect = "coordinates"),
               "constant within")
  expect_error(fit_model(d, responses = "biomass", spatial_effect = "coordinates",
                         spatial_options = list(k = 2)), "k")
  expect_error(fit_model(d, responses = "biomass", spatial_effect = "coordinates",
                         spatial_options = list(foo = 1)), "spatial_options")
  d3 <- d
  d3$latitude <- as.character(d3$latitude)
  expect_error(fit_model(d3, responses = "biomass", spatial_effect = "coordinates"),
               "numeric")
})

test_that("spatial models compile and run, alone and with a locality effect", {
  skip_if_not_installed("nimble")
  d <- sp_data("autocorrelation")
  f1 <- bm_test_fit(d, responses = "biomass", spatial_effect = "autocorrelation",
                    random_effects = "year")
  expect_true(all(c("sd_sp_biomass", "phi_idx_sp_biomass") %in% f1$diagnostics$parameter))
  f2 <- bm_test_fit(d, responses = c("biomass", "richness"),
                    spatial_effect = "coordinates",
                    random_effects = c("year", "locality"))
  expect_true(all(c("sd_sp_biomass", "sd_sp_richness", "sd_locality_biomass") %in%
                    f2$diagnostics$parameter))
  expect_false(any(grepl("^(w_sp|v_sp)", summary(f2)$table$parameter)))
  # coordinates passed explicitly
  d$lon2 <- d$longitude
  d$lat2 <- d$latitude
  f3 <- bm_test_fit(d, responses = "biomass", spatial_effect = "coordinates",
                    random_effects = "year", coordinates = c("lon2", "lat2"))
  expect_s3_class(f3, "bm_fit")
})

test_that("spatial predictions at new sites are finite and deterministic when requested", {
  skip_if_not_installed("nimble")
  d <- sp_data("autocorrelation")
  f <- bm_test_fit(d, responses = "biomass", spatial_effect = "autocorrelation",
                   random_effects = "year")
  new_xy <- cbind(c(8, 12), c(60, 65))
  a <- bm_spatial_draws(f, "biomass", new_xy, simulate = FALSE)
  b <- bm_spatial_draws(f, "biomass", new_xy, simulate = FALSE)
  expect_equal(dim(a), c(nrow(as.matrix(f$samples)), 2L))
  expect_true(all(is.finite(a)))
  expect_equal(a, b)
  s <- bm_spatial_draws(f, "biomass", new_xy, simulate = TRUE)
  expect_gt(mean(apply(s, 2, stats::sd)), mean(apply(a, 2, stats::sd)))
})

test_that("NNGP recovers the field and predicts held-out sites (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(
    n_year = 4, n_locality = 70, n_rep = 3, seed = 55, spatial = "autocorrelation",
    spatial_sd = c(biomass = 8, richness = 100), spatial_range_km = 250,
    sd_locality = c(biomass = 0.01, richness = 1),
    sd_replication = c(biomass = 0.01, richness = 1)
  )
  truth <- attr(d, "truth")
  locs <- unique(d$locality)
  held <- locs[56:70]
  fit_d <- d[!d$locality %in% held, ]
  fit <- bm_test_fit(fit_d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     spatial_effect = "autocorrelation", random_effects = "year",
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  est <- fit$diagnostics
  row <- est[est$parameter == "sd_sp_biomass", ]
  expect_true(row$q2.5 <= 8 && 8 <= row$q97.5)
  row <- est[est$parameter == "sigma_biomass", ]
  expect_true(row$q2.5 <= truth$sigma[["biomass"]] && truth$sigma[["biomass"]] <= row$q97.5)
  sites <- fit$prep$spatial$site_names
  w <- colMeans(as.matrix(fit$samples)[, paste0("w_sp_biomass[", seq_along(sites), "]")])
  expect_gt(cor(w, truth$spatial_field$biomass[sites]), 0.9)
  # Held-out sites
  hxy <- unique(d[d$locality %in% held, c("locality", "longitude", "latitude")])
  pred <- colMeans(bm_spatial_draws(fit, "biomass", hxy[, 2:3], simulate = FALSE))
  expect_gt(cor(pred, truth$spatial_field$biomass[hxy$locality]), 0.6)
})

test_that("coordinate spline recovers a smooth surface (long run)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  d <- simulate_bm_data(
    n_year = 4, n_locality = 70, n_rep = 3, seed = 77, spatial = "coordinates",
    spatial_sd = c(biomass = 6, richness = 60),
    sd_locality = c(biomass = 0.01, richness = 1),
    sd_replication = c(biomass = 0.01, richness = 1)
  )
  truth <- attr(d, "truth")
  fit <- bm_test_fit(d, responses = "biomass", fixed_effects = c("x1", "x2"),
                     spatial_effect = "coordinates", random_effects = "year",
                     n_chains = 3, n_iter = 20000, n_burnin = 10000, n_thin = 10)
  sites <- fit$prep$spatial$site_names
  est <- colMeans(bm_spatial_draws(fit, "biomass", fit$prep$spatial$xy))
  tf <- truth$spatial_field$biomass[sites]
  expect_gt(cor(est, tf), 0.9)
  expect_lt(sqrt(mean(((est - mean(est)) - (tf - mean(tf)))^2)), 2)
})
