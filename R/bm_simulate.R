#' Simulate a biomass and richness data set with known parameters
#'
#' Generates a balanced data set from the null model (plus optional covariate
#' effects) using the truncated normal likelihood. Used for examples and for
#' validating that parameters are recovered. The true values are returned in
#' `attr(, "truth")`.
#'
#' @param n_year,n_locality,n_rep Number of years, localities and replicates
#'   per locality-year.
#' @param beta0 Named numeric, intercepts of the untruncated normal for
#'   `biomass` and `richness`.
#' @param beta Named list, covariate effects (per sd of `x1`, `x2`) by response.
#' @param sd_year,sd_locality,sd_replication,sigma Named numeric vectors by
#'   response.
#' @param trend Named numeric, linear trend per year by response (0 = none).
#' @param cyclic_amplitude Named numeric, amplitude of a seasonal
#'   `amplitude * sin(2 * pi * day / 365)` effect by response (0 = none).
#'   Day of year is stored in `start_date_julian`.
#' @param spatial `"none"`, `"autocorrelation"` (exponential Gaussian process
#'   over localities) or `"coordinates"` (a smooth function of the planar
#'   coordinates). The true site-level fields are stored in
#'   `attr(, "truth")$spatial_field`, centred to mean zero.
#' @param spatial_sd Named numeric, spatial marginal sd (or amplitude) by
#'   response.
#' @param spatial_range_km Length scale `ell` (km) of the exponential
#'   correlation `exp(-d / ell)` for the autocorrelated field.
#' @param interaction Logical. If `TRUE` (needs a `spatial` type), adds a spatial
#'   field that modulates the temporal term: a spatially varying slope
#'   (`trend + field`) per year when `trend` is non-zero for a response, otherwise a
#'   spatially varying amplitude (`(1 + field) * cyclic effect`) of the
#'   seasonal effect.
#' @param interaction_sd Named numeric, sd (or amplitude) of the interaction
#'   field by response: slope units per year, or a dimensionless multiplier.
#'   Required when `interaction = TRUE`.
#' @param loadings Named numeric. For each named component, the richness effect
#'   is `loading * biomass effect` (a shared latent process with loading):
#'   components `"x1"`, `"x2"` (coefficients), `"year"`, `"locality"`,
#'   `"replication"` (random effects; the richness sd of that effect is then
#'   `abs(loading) * biomass sd`), `"trend"`, `"cyclic"`, `"spatial"` and
#'   `"interaction"` (fields). Unnamed components are generated independently.
#' @param rho Named numeric by random effect (`"year"`, `"locality"`,
#'   `"replication"`): correlation between the biomass and richness random
#'   effects (each keeping its own sd).
#' @param family Distribution of each response, `"truncnorm"` (default) or
#'   `"negbinomial"` for richness. For a negative binomial response the linear
#'   predictor is on the **log scale**: its `beta0`, `beta`, `sd_*`, `trend`,
#'   `cyclic_amplitude` and `spatial_sd` (and the loadings applied to them) are
#'   log-scale values, e.g. `beta0 = log(500)`, and `sigma` is not used.
#' @param nb_size Size (dispersion) `phi` of the negative binomial, with variance
#'   `mean + mean^2 / phi`.
#' @param seed Integer seed.
#' @return A data frame with the default columns of [bm_default_columns()]
#'   plus covariates `x1`, `x2`, with attribute `truth`.
#' @export
#' @examples
#' sim <- simulate_bm_data(n_year = 3, n_locality = 5, n_rep = 2, seed = 7)
#' head(sim)
#' attr(sim, "truth")$beta0
#' @md
simulate_bm_data <- function(n_year = 5, n_locality = 30, n_rep = 4,
                             beta0 = c(biomass = 20, richness = 500),
                             beta = list(biomass = c(x1 = 5, x2 = 0),
                                         richness = c(x1 = 0, x2 = 80)),
                             sd_year = c(biomass = 3, richness = 40),
                             sd_locality = c(biomass = 6, richness = 100),
                             sd_replication = c(biomass = 4, richness = 50),
                             sigma = c(biomass = 5, richness = 80),
                             trend = c(biomass = 0, richness = 0),
                             cyclic_amplitude = c(biomass = 0, richness = 0),
                             spatial = c("none", "autocorrelation", "coordinates"),
                             spatial_sd = c(biomass = 6, richness = 100),
                             spatial_range_km = 300,
                             interaction = FALSE,
                             interaction_sd = NULL,
                             family = c(biomass = "truncnorm", richness = "truncnorm"),
                             nb_size = 5,
                             loadings = NULL,
                             rho = NULL,
                             seed = NULL) {
  spatial <- match.arg(spatial)
  sim <- function() {
    d <- expand.grid(rep = seq_len(n_rep), locality = seq_len(n_locality),
                     year = seq_len(n_year))
    d$year <- 2019L + d$year
    d$locality <- sprintf("loc_%03d", d$locality)
    d$sampling_name <- sprintf("%s_%d_sampling-%02d", d$locality, d$year, d$rep)
    loc <- unique(d$locality)
    xy <- data.frame(locality = loc, longitude = stats::runif(length(loc), 5, 15),
                     latitude = stats::runif(length(loc), 58, 68))
    d <- merge(d, xy, by = "locality", sort = FALSE)
    field <- list(biomass = 0, richness = 0)
    field_int <- list(biomass = 0, richness = 0)
    if (spatial != "none") {
      xyk <- bm_plane(as.matrix(xy[, c("longitude", "latitude")]),
                      list(lonlat = TRUE, x0 = 10, y0 = 63, scale = 1))
      Dkm <- bm_distance_km(as.matrix(xy[, 2:3]), as.matrix(xy[, 2:3]), TRUE)
      gen_field <- function(sdv, range_km) {
        if (spatial == "autocorrelation") {
          C <- sdv^2 * exp(-Dkm / range_km) + diag(1e-8, nrow(xy))
          v <- as.numeric(t(chol(C)) %*% stats::rnorm(nrow(xy)))
        } else {
          ph <- stats::runif(2, 0, 2 * pi)
          v <- sdv * (sin(xyk[, 1] / 400 + ph[1]) + cos(xyk[, 2] / 400 + ph[2]))
        }
        stats::setNames(v - mean(v), xy$locality)
      }
      for (r in names(field)) field[[r]] <- gen_field(spatial_sd[[r]], spatial_range_km)
      if (interaction) {
        for (r in names(field_int)) {
          field_int[[r]] <- gen_field(interaction_sd[[r]], spatial_range_km)
        }
      }
    }
    d$start_date_julian <- as.integer(pmin(364, pmax(1, round(
      120 + (d$rep - 1) * 150 / max(n_rep - 1, 1) +
        stats::runif(nrow(d), -12, 12)))))
    d$x1 <- stats::rnorm(nrow(d))
    d$x2 <- stats::rnorm(nrow(d))
    re <- function(ids, s) stats::rnorm(length(unique(ids)), 0, s)[as.integer(factor(ids))]
    truth <- list(beta0 = beta0, beta = beta, sd_year = sd_year,
                  sd_locality = sd_locality, sd_replication = sd_replication,
                  sigma = sigma, family = family, nb_size = nb_size, spatial = spatial, spatial_field = field, interaction_field = field_int,
                  interaction = interaction, loadings = loadings, rho = rho,
                  spatial_sd = spatial_sd, spatial_range_km = spatial_range_km,
                  trend = trend, cyclic_amplitude = cyclic_amplitude,
                  cyclic_fun = function(doy, r) {
                    cyclic_amplitude[[r]] * sin(2 * pi * doy / 365)
                  })
    u_b <- list()
    ld <- function(nm) if (!is.null(loadings) && !is.na(loadings[nm])) loadings[[nm]]
    for (r in c("biomass", "richness")) {
      sh <- r == "richness"
      # Effects of richness that are shared are `loading * biomass effect`
      scal <- function(nm, value_b, value_r) {
        l <- if (sh) ld(nm)
        if (is.null(l)) value_r else l * value_b
      }
      b1 <- scal("x1", beta$biomass[["x1"]], beta[[r]][["x1"]])
      b2 <- scal("x2", beta$biomass[["x2"]], beta[[r]][["x2"]])
      tr <- scal("trend", trend[["biomass"]], trend[[r]])
      amp <- scal("cyclic", cyclic_amplitude[["biomass"]], cyclic_amplitude[[r]])
      fld <- scal("spatial", field$biomass, field[[r]])
      fint <- scal("interaction", field_int$biomass, field_int[[r]])
      # Record the effective (possibly shared) values as the truth
      truth$beta[[r]][["x1"]] <- b1
      truth$beta[[r]][["x2"]] <- b2
      truth$trend[[r]] <- tr
      truth$cyclic_amplitude[[r]] <- amp
      truth$spatial_field[[r]] <- fld
      truth$interaction_field[[r]] <- fint
      u <- list()
      for (e in c("year", "locality", "replication")) {
        ids <- switch(e, year = d$year, locality = d$locality, replication = d$rep)
        sdv <- list(year = sd_year, locality = sd_locality,
                    replication = sd_replication)[[e]][[r]]
        eps <- re(ids, sdv)
        if (!sh) {
          u[[e]] <- eps
          u_b[[e]] <- eps
        } else if (!is.null(ld(e))) {
          u[[e]] <- ld(e) * u_b[[e]]
        } else if (!is.null(rho) && !is.na(rho[e])) {
          sdb <- list(year = sd_year, locality = sd_locality,
                      replication = sd_replication)[[e]][["biomass"]]
          u[[e]] <- rho[[e]] * sdv / sdb * u_b[[e]] + sqrt(1 - rho[[e]]^2) * eps
        } else {
          u[[e]] <- eps
        }
      }
      # Realised (not nominal) site effects, for validating correlations
      truth$locality_effect[[r]] <- tapply(u$locality, d$locality, `[`, 1)
      cyc <- amp * sin(2 * pi * d$start_date_julian / 365)
      mu <- beta0[[r]] + b1 * d$x1 + b2 * d$x2 +
        (tr + (if (interaction) fint[d$locality] else 0)) * (d$year - mean(d$year)) +
        (1 + (if (interaction && tr == 0) fint[d$locality] else 0)) * cyc +
        (if (spatial == "none") 0 else fld[d$locality]) +
        u$year + u$locality + u$replication
      if (identical(unname(family[[r]]), "negbinomial")) {
        y <- stats::rnbinom(nrow(d), size = nb_size, mu = exp(mu))
      } else {
        lo <- stats::pnorm(0, mu, sigma[[r]])
        y <- stats::qnorm(stats::runif(nrow(d), lo, 1), mu, sigma[[r]])
        y <- pmax(y, 0)
      }
      d[[if (r == "biomass") "sum_wet_weight" else "no_species"]] <- y
    }
    d$rep <- NULL
    d <- d[order(d$year, d$locality, d$sampling_name), ]
    rownames(d) <- NULL
    attr(d, "truth") <- truth
    d
  }
  if (is.null(seed)) sim() else bm_with_seed(seed, sim())
}
