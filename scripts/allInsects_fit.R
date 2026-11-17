# Joint biomass (truncated normal) and richness (negative binomial) model for
# all insects: NNGP spatial field + cyclic seasonal spline, both shared.
# Requires `all_insect_data` (from 00_masterScript.R data preparation).

library(BINTS)
library(dplyr)
library(ggplot2)

out_dir <- "scripts/allInsects_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- all_insect_data |>
  filter(!is.na(sum_wet_weight), sum_wet_weight >= 0, !is.na(no_species),
         !is.na(mean_temp), !is.na(mean_precip), !is.na(start_date_julian),
         !is.na(longitude), !is.na(latitude)) |>
  mutate(no_species = round(no_species)) |>
  as.data.frame()

fit <- fit_model(
  dat,
  responses = c("biomass", "richness"),
  family = c(biomass = "truncnorm", richness = "negbinomial"),
  fixed_effects = c("mean_temp", "mean_precip"),
  random_effects = c("year", "locality", "replication"),
  temporal_effect = "cyclic",
  spatial_effect = "autocorrelation",
  spatial_options = list(m = 15),
  shared_spatial_effect = TRUE,
  shared_temporal_effect = TRUE,
  n_chains = get0("n_chains_fit", ifnotfound = 4),
  n_iter = get0("n_iter_fit", ifnotfound = 6000),
  n_burnin = get0("n_burnin_fit", ifnotfound = 3000),
  n_thin = get0("n_thin_fit", ifnotfound = 3),
  seed = 4817, parallel = get0("parallel_fit", ifnotfound = TRUE)
)
saveRDS(fit, file.path(out_dir, "fit.rds"))

ana <- function(fit, n_draws = 500, cell_km = 25) {
  # Posterior predictive checks
  ppc <- posterior_predictive_check(fit, n_draws = n_draws, seed = 2951)
  saveRDS(ppc, file.path(out_dir, "ppc.rds"))
  for (r in c("biomass", "richness")) {
    ggsave(file.path(out_dir, paste0("ppc_", r, ".png")),
           plot(ppc, response = r), width = 8, height = 5)
  }

  # Residuals and observed vs predicted
  res <- residuals(fit)
  saveRDS(res, file.path(out_dir, "residuals.rds"))
  for (r in c("biomass", "richness")) {
    ggsave(file.path(out_dir, paste0("obs_pred_", r, ".png")),
           plot_observed_predicted(fit, response = r), width = 6, height = 6)
    ggsave(file.path(out_dir, paste0("residuals_", r, ".png")),
           plot_residuals(fit, covariates = c("mean_temp", "mean_precip", "start_date_julian"),
                          response = r), width = 9, height = 6)
  }

  # Prediction over the region: hull of sampled localities, covariates held at
  # their means (no covariate surfaces for the grid), mid-season, last year
  pts <- sf::st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(25833)
  hull <- sf::st_buffer(sf::st_convex_hull(sf::st_union(pts)), 20000)
  grid <- sf::st_make_grid(hull, cellsize = cell_km * 1000, what = "centers")
  grid <- grid[sf::st_intersects(grid, hull, sparse = FALSE)[, 1]] |>
    sf::st_transform(4326)
  xy <- sf::st_coordinates(grid)
  gdf <- data.frame(longitude = xy[, 1], latitude = xy[, 2], year = max(dat$year),
                    start_date_julian = 200, mean_temp = mean(dat$mean_temp),
                    mean_precip = mean(dat$mean_precip))
  gsf <- sf::st_as_sf(gdf, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  pred <- predict(fit, gsf, random = c(year = "zero", locality = "zero", replication = "zero"),
                  n_draws = 200)
  saveRDS(pred, file.path(out_dir, "pred_grid.rds"))
  for (r in c("biomass", "richness")) {
    ggsave(file.path(out_dir, paste0("map_", r, ".png")),
           plot_spatial_predictions(fit, pred, r, uncertainty = "width"),
           width = 7, height = 8)
  }
  invisible(list(ppc = ppc, res = res, pred = pred))
}

out <- ana(fit)
