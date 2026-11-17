# =============================================================================
# create_spatial_shared_model
# Spatial-only model with a shared latent spatial process
# Standalone model builder; not a wrapper around another model function.
# =============================================================================

#' Create a NIMBLE Model with Shared Fixed, Random, and Spatial Effects
#' for Insect Biomass and Diversity
#'
#' Builds a NIMBLE model for insect biomass, diversity, or both responses.
#' Continuous covariate coefficients, categorical covariate effects, enabled
#' year/locality/replicate random effects, and spatial effects are shared
#' between biomass and diversity. Response-specific intercepts, observation
#' distributions, and distribution-specific parameters remain separate.
#'
#' @param modelCovs Character vector of covariate names to include.
#' @param response Response to model: `"biomass"`, `"diversity"`, or `"both"`.
#'   Defaults to `"biomass"`.
#' @param biomass_distribution Biomass distribution: `"normal"` or
#'   `"lognormal"`. Defaults to `"normal"`.
#' @param diversity_distribution Diversity distribution: `"poisson"` or
#'   `"negative_binomial"`. Defaults to `"poisson"`.
#' @param covariate_types Optional named character vector classifying each
#'   covariate as `"continuous"` or `"categorical"`. If omitted, all covariates
#'   are treated as continuous. Names must match `modelCovs`.
#' @param random_year Logical; include a shared year random effect if `TRUE`.
#'   Defaults to `FALSE`.
#' @param random_locality Logical; include a shared locality random effect if
#'   `TRUE`. Defaults to `FALSE`.
#' @param random_replicate Logical; include a shared replicate random effect if
#'   `TRUE`. Defaults to `FALSE`.
#' @param temporal_effect Temporal effect: `"none"`, `"linear"`, or
#'   `"cubic_cyclic"`. Defaults to `"none"`. Temporal effects remain
#'   response-specific unless `temporal_shared = TRUE` for a cubic cyclic
#'   spline.
#' @param temporal_basis Numeric matrix of cyclic spline basis values. Required
#'   when `temporal_effect = "cubic_cyclic"`; rows must align with observations.
#' @param temporal_penalty Numeric square penalty matrix for the cyclic spline.
#'   Required for `temporal_effect = "cubic_cyclic"` and must have dimensions
#'   equal to `ncol(temporal_basis)`.
#' @param temporal_penalty_tol Positive finite threshold for retaining positive
#'   penalty eigenvalues. Defaults to `1e-8`.
#' @param temporal_shared Logical; use a shared latent temporal spline process
#'   if `TRUE`. Currently supported only for `"cubic_cyclic"`. Defaults to
#'   `FALSE`.
#' @param spatial_effect Spatial effect: `"none"`, `"coordinates"`, or
#'   `"autocorrelation"`. Defaults to `"autocorrelation"`.
#' @param spatial_shared Must be `TRUE`. The same spatial effect is used in
#'   both response predictors when both responses are included. Defaults to
#'   `TRUE`; `FALSE` is rejected.
#' @param spatial_range_prior Positive numeric vector of length two giving the
#'   lower and upper bounds for the spatial range prior. The upper bound must
#'   exceed the lower bound. Defaults to `c(1000, 500000)`.
#' @param spatial_sigma_prior Positive finite scalar controlling the prior
#'   scale for the spatial standard deviation. Defaults to `1`.
#' @param output_dir Directory for the generated model file. Defaults to
#'   `temp_dir`, which must exist in the calling environment.
#' @param model_name Base name for the generated R file. Defaults to
#'   `"create_spatial_shared_model"`; the `.R` extension is added automatically.
#'
#' @details
#' Continuous covariates use a single coefficient vector shared by the
#' included responses. Each categorical covariate has one shared vector of
#' level effects and one precision parameter. Enabled year, locality, and
#' replicate effects likewise use one shared effect vector and precision
#' parameter per grouping factor.
#'
#' With `spatial_effect = "coordinates"`, the longitude and latitude
#' coefficients (including their quadratic terms) are shared. With
#' `spatial_effect = "autocorrelation"`, the responses use the same latent
#' spatial Gaussian process with exponential covariance
#' \deqn{\Sigma_{ij} = \sigma^2 \exp(-d_{ij}/\rho),}
#' where `d_ij` is the distance between locations, `rho` is the spatial range,
#' and `sigma` is the spatial standard deviation. The shared spatial effect is
#' added directly to each included response predictor; no response-specific
#' spatial loading is estimated.
#'
#' Biomass and diversity retain separate intercepts and observation models.
#' Biomass uses a normal or lognormal likelihood; diversity uses a Poisson or
#' negative-binomial likelihood. Temporal terms are not shared by default:
#' `temporal_shared = TRUE` enables sharing for the cubic cyclic spline only.
#'
#' The generated NIMBLE model expects required data, dimensions, indices, and
#' constants to be supplied when the model is instantiated. These may include
#' `X`, `longitude`, `latitude`, `distance_matrix`, `spatial_zero`,
#' `spatial_id`, and the enabled grouping indices `year_id`, `locality_id`,
#' and `replicate_id`. For cyclic splines, use the returned
#' `temporal_constants` list as NIMBLE constants.
#'
#' @return Invisibly returns a list containing the generated model text and
#'   output file path, response and distribution settings, covariate
#'   classifications, random-effect flags, temporal spline information, and
#'   spatial configuration and prior settings.
#'
#' @examples
#' \dontrun{
#' model_info <- create_insect_spatial_shared_model(
#'   modelCovs = c("temperature", "vegetation"),
#'   response = "both",
#'   covariate_types = c(
#'     temperature = "continuous",
#'     vegetation = "categorical"
#'   ),
#'   random_year = TRUE,
#'   random_locality = TRUE,
#'   spatial_effect = "autocorrelation",
#'   spatial_shared = TRUE,
#'   output_dir = tempdir(),
#'   model_name = "insect_shared_spatial"
#' )
#' }
#'
#' @seealso [nimble::nimbleCode()]
#' @export

create_insect_spatial_shared_model <- function(

  modelCovs,

  response = "biomass",

  # --------------------------------------------------------------------------
  # Response distributions
  # --------------------------------------------------------------------------

  biomass_distribution = "normal",

  diversity_distribution = "poisson",

  # --------------------------------------------------------------------------
  # Covariate specification
  # --------------------------------------------------------------------------

  covariate_types = NULL,
  share_fixed_effects = FALSE,
  share_hyperparameters = FALSE,
  
  # --------------------------------------------------------------------------
  # Random effects
  # --------------------------------------------------------------------------

  random_year = FALSE,

  random_locality = FALSE,

  random_replicate = FALSE,

  # --------------------------------------------------------------------------
  # Temporal effect
  # --------------------------------------------------------------------------

  temporal_effect = "none",

  temporal_basis = NULL,

  temporal_penalty = NULL,

  temporal_penalty_tol = 1e-8,

  # If TRUE, use one cyclic-spline process for both responses.
  # For a two-response model, the biomass loading is fixed to 1 and the
  # diversity loading is estimated. Shared temporal effects are currently
  # implemented for cubic_cyclic splines.
  temporal_shared = FALSE,

  # --------------------------------------------------------------------------
  # Spatial effect
  # --------------------------------------------------------------------------

  spatial_effect = c("autocorrelation", "coordinates"),

  spatial_shared = TRUE,

  # --------------------------------------------------------------------------
  # Spatial priors
  # --------------------------------------------------------------------------

  spatial_range_prior = c(
    1000,
    500000
  ),

  spatial_sigma_prior = 1,

  # --------------------------------------------------------------------------
  # Output
  # --------------------------------------------------------------------------

  output_dir = temp_dir,

  model_name = "create_spatial_shared_model"

) {

  # Fixed structure for this standalone model builder: Spatial-only model with a shared latent spatial process.
  effect_structure <- "spatial"


  # ============================================================================
  # 1. Validate response
  # ============================================================================

  response <- match.arg(
    response,
    choices = c(
      "biomass",
      "diversity",
      "both"
    )
  )


  # ============================================================================
  # 2. Validate response distributions
  # ============================================================================

  biomass_distribution <- match.arg(
    biomass_distribution,
    choices = c(
      "normal",
      "lognormal"
    )
  )

  diversity_distribution <- match.arg(
    diversity_distribution,
    choices = c(
      "poisson",
      "negative_binomial"
    )
  )


  # ============================================================================
  # 3. Validate temporal effect
  # ============================================================================

  temporal_effect <- match.arg(
    temporal_effect
  )

  spatial_effect <- match.arg(spatial_effect)

  if (!is.logical(temporal_shared) || length(temporal_shared) != 1L ||
      is.na(temporal_shared)) {
    stop("'temporal_shared' must be TRUE or FALSE.", call. = FALSE)
  }

  allowed_structures <- c("none", "temporal", "spatial", "spatiotemporal")

  if (is.null(effect_structure)) {
    has_temporal <- temporal_effect != "none"
    has_spatial <- spatial_effect != "none"
    effect_structure <- if (has_temporal && has_spatial) {
      "spatiotemporal"
    } else if (has_temporal) {
      "temporal"
    } else if (has_spatial) {
      "spatial"
    } else {
      "none"
    }
  } else {
    effect_structure <- match.arg(effect_structure, allowed_structures)
  }

  # Make effect_structure authoritative. This prevents a temporal or spatial
  # term from being included accidentally in a candidate model.
  if (effect_structure %in% c("none", "spatial")) {
    temporal_effect <- "none"
  } else if (temporal_effect == "none") {
    stop(
      "For effect_structure = '", effect_structure,
      "', choose temporal_effect = 'linear' or 'cubic_cyclic'.",
      call. = FALSE
    )
  }

  if (effect_structure %in% c("none", "temporal")) {
    spatial_effect <- "none"
  } else if (spatial_effect == "none") {
    stop(
      "For effect_structure = '", effect_structure,
      "', choose spatial_effect = 'coordinates' or 'autocorrelation'.",
      call. = FALSE
    )
  }


    temporal_basis_penalized <- NULL

    temporal_penalty_eigenvalues <- NULL

    temporal_penalty_eigenvectors <- NULL

    n_temporal_basis <- 0L

  


  # ============================================================================
  # 4. Validate spatial effect
  # ============================================================================

  spatial_effect <- match.arg(
    spatial_effect
  )


  if (
    !is.logical(spatial_shared) ||
    length(spatial_shared) != 1
  ) {

    stop(
      "'spatial_shared' must be TRUE or FALSE.",
      call. = FALSE
    )

  }


  if (
    spatial_effect != "autocorrelation" &&
    !isTRUE(spatial_shared)
  ) {

    warning(
      "'spatial_shared' is only used for spatial autocorrelation. ",
      "It will be ignored for spatial_effect = '",
      spatial_effect,
      "'.",
      call. = FALSE
    )

  }


  # ============================================================================
  if (!isTRUE(spatial_shared)) {
    stop(
      "'spatial_shared' must be TRUE for create_insect_spatial_shared_model().",
      call. = FALSE
    )
  }

  # ===========================================================================
  # 5. Validate spatial priors
  # ============================================================================

  if (
    !is.numeric(spatial_range_prior) ||
    length(spatial_range_prior) != 2 ||
    any(!is.finite(spatial_range_prior)) ||
    spatial_range_prior[1] <= 0 ||
    spatial_range_prior[2] <= spatial_range_prior[1]
  ) {

    stop(
      "'spatial_range_prior' must contain two positive values ",
      "with upper > lower.",
      call. = FALSE
    )

  }


  if (
    !is.numeric(spatial_sigma_prior) ||
    length(spatial_sigma_prior) != 1 ||
    !is.finite(spatial_sigma_prior) ||
    spatial_sigma_prior <= 0
  ) {

    stop(
      "'spatial_sigma_prior' must be a single positive finite number.",
      call. = FALSE
    )

  }


  # ============================================================================
  # 6. Determine included responses
  # ============================================================================

  include_biomass <- response %in% c(
    "biomass",
    "both"
  )

  include_diversity <- response %in% c(
    "diversity",
    "both"
  )


  # ============================================================================
  # 7. Validate covariate types
  # ============================================================================

  if (length(modelCovs) > 0) {

    if (is.null(covariate_types)) {

      covariate_types <- setNames(
        rep(
          "continuous",
          length(modelCovs)
        ),
        modelCovs
      )

    } else {

      invalid_types <- setdiff(
        unique(covariate_types),
        c(
          "continuous",
          "categorical"
        )
      )

      if (length(invalid_types) > 0) {

        stop(
          "covariate_types must contain only ",
          "'continuous' or 'categorical'.",
          call. = FALSE
        )

      }


      missing_types <- setdiff(
        modelCovs,
        names(covariate_types)
      )

      if (length(missing_types) > 0) {

        stop(
          "No covariate type supplied for: ",
          paste(
            missing_types,
            collapse = ", "
          ),
          call. = FALSE
        )

      }

    }

  } else {

    covariate_types <- character(0)

  }


  # ============================================================================
  # 8. Separate continuous and categorical covariates
  # ============================================================================

  continuous_covariates <- modelCovs[
    covariate_types[modelCovs] == "continuous"
  ]

  categorical_covariates <- modelCovs[
    covariate_types[modelCovs] == "categorical"
  ]


  # ============================================================================
  # 9. Shared fixed effects
  # ============================================================================
 # ============================================================================
  # 9. Fixed effects
  # ============================================================================

  n_continuous <- length(
    continuous_covariates
  )


  if (n_continuous > 0) {

    if(share_fixed_effects){

    fixed_effect_all <- paste0(
      "beta[",
      seq_along(continuous_covariates),
      "] * X[i, ",
      seq_along(continuous_covariates),
      "]",
      collapse = " + "
    )
      
    fixed_effect_biomass <- fixed_effect_all
    fixed_effect_diversity <- fixed_effect_all
      
    } else {

    fixed_effect_biomass <- paste0(
      "beta_biomass[",
      seq_along(continuous_covariates),
      "] * X[i, ",
      seq_along(continuous_covariates),
      "]",
      collapse = " + "
    )

    fixed_effect_diversity <- paste0(
      "beta_diversity[",
      seq_along(continuous_covariates),
      "] * X[i, ",
      seq_along(continuous_covariates),
      "]",
      collapse = " + "
    )
  }

    } else {

    fixed_effect_biomass <- NULL

    fixed_effect_diversity <- NULL

  }


  # 10. Initialise random effects
  # ============================================================================

  random_effect_code <- character(0)

  biomass_random_terms <- character(0)

  diversity_random_terms <- character(0)


  # ============================================================================
  # 11. Initialise spatial effects
  # ============================================================================

  spatial_code <- character(0)

  biomass_spatial_terms <- character(0)

  diversity_spatial_terms <- character(0)


  # ============================================================================
  # 12. Spatial coordinate effects
  # ============================================================================

  if (spatial_effect == "coordinates") {

    spatial_code <- c(
      spatial_code,

      '
# =============================================================================
# Shared spatial coordinate effect
# =============================================================================

beta_longitude_spatial ~ dnorm(
  0,
  0.001
)

beta_latitude_spatial ~ dnorm(
  0,
  0.001
)

beta_longitude_spatial_squared ~ dnorm(
  0,
  0.001
)

beta_latitude_spatial_squared ~ dnorm(
  0,
  0.001
)
'
    )


if (include_biomass) {

  biomass_spatial_terms <- c(
    biomass_spatial_terms,

    "beta_longitude_spatial * longitude[i]",

    "beta_longitude_spatial_squared * pow(longitude[i], 2)",

    "beta_latitude_spatial * latitude[i]",

    "beta_latitude_spatial_squared * pow(latitude[i], 2)"
  )

}


if (include_diversity) {

  diversity_spatial_terms <- c(
    diversity_spatial_terms,

    "beta_longitude_spatial * longitude[i]",

    "beta_longitude_spatial_squared * pow(longitude[i], 2)",

    "beta_latitude_spatial * latitude[i]",

    "beta_latitude_spatial_squared * pow(latitude[i], 2)"
  )

}

  }


  # ============================================================================
  # 13. Spatial autocorrelation
  # ============================================================================

  if (spatial_effect == "autocorrelation") {

    # --------------------------------------------------------------------------
    # Shared spatial process
    # --------------------------------------------------------------------------

    if (spatial_shared) {

      spatial_code <- c(
        spatial_code,

        '
# =============================================================================
# Shared latent spatial autocorrelation
# =============================================================================

rho_spatial ~ dunif(
  spatial_range_lower,
  spatial_range_upper
)

sigma_spatial ~ T(
  dnorm(
    0,
    spatial_sigma_precision
  ),
  0,
  Inf
)

for (i in 1:n_spatial) {

  for (j in 1:n_spatial) {

    spatial_covariance[i, j] <-
      sigma_spatial^2 *
      exp(
        -distance_matrix[i, j] /
        rho_spatial
      )

  }

}

spatial_effect[1:n_spatial] ~ dmnorm(
  mean = spatial_zero[1:n_spatial],
  cov = spatial_covariance[
    1:n_spatial,
    1:n_spatial
  ]
)

# The same latent spatial effect is added to each included response predictor.
'
      )



if (include_biomass) {
  biomass_spatial_terms <- c(
    biomass_spatial_terms,
    "spatial_effect[spatial_id[i]]"
  )
}

if (include_diversity) {
  diversity_spatial_terms <- c(
    diversity_spatial_terms,
    "spatial_effect[spatial_id[i]]"
  )
}

    }


  # --------------------------------------------------------------------------
  # Separate spatial processes
  # --------------------------------------------------------------------------

  if (!spatial_shared) {

    if (include_biomass) {

      spatial_code <- c(
        spatial_code,

        '
# =============================================================================
# Biomass latent spatial autocorrelation
# =============================================================================

rho_spatial_biomass ~ dunif(
  spatial_range_lower,
  spatial_range_upper
)

sigma_spatial_biomass ~ T(
  dnorm(
    0,
    spatial_sigma_precision
  ),
  0,
  Inf
)

for (i in 1:n_spatial) {

  for (j in 1:n_spatial) {

    spatial_covariance_biomass[i, j] <-
      sigma_spatial_biomass^2 *
      exp(
        -distance_matrix[i, j] /
        rho_spatial_biomass
      )

  }

}

spatial_effect_biomass[1:n_spatial] ~ dmnorm(
  mean = spatial_zero_biomass[1:n_spatial],
  cov = spatial_covariance_biomass[
    1:n_spatial,
    1:n_spatial
  ]
)

lambda_spatial_biomass ~ dnorm(
  0,
  0.001
)
'
      )


biomass_spatial_terms <- c(
  biomass_spatial_terms,

  "lambda_spatial_biomass * spatial_effect_biomass[spatial_id[i]]"
)

    }


  if (include_diversity) {

    spatial_code <- c(
      spatial_code,

      '
# =============================================================================
# Diversity latent spatial autocorrelation
# =============================================================================

rho_spatial_diversity ~ dunif(
  spatial_range_lower,
  spatial_range_upper
)

sigma_spatial_diversity ~ T(
  dnorm(
    0,
    spatial_sigma_precision
  ),
  0,
  Inf
)

for (i in 1:n_spatial) {

  for (j in 1:n_spatial) {

    spatial_covariance_diversity[i, j] <-
      sigma_spatial_diversity^2 *
      exp(
        -distance_matrix[i, j] /
        rho_spatial_diversity
      )

  }

}

spatial_effect_diversity[1:n_spatial] ~ dmnorm(
  mean = spatial_zero_diversity[1:n_spatial],
  cov = spatial_covariance_diversity[
    1:n_spatial,
    1:n_spatial
  ]
)

lambda_spatial_diversity ~ dnorm(
  0,
  0.001
)
'
    )


diversity_spatial_terms <- c(
  diversity_spatial_terms,

  "lambda_spatial_diversity * spatial_effect_diversity[spatial_id[i]]"
)

  }

  }

  }


# ============================================================================
# 14. Shared categorical covariate effects
# ============================================================================

if (length(categorical_covariates) > 0) {

  for (cov in categorical_covariates) {

    cov_clean <- make.names(cov)


    if(share_hyperparameters){
        paste0(
          '
# =============================================================================
# Hyperparameters for: ', cov, '
# =============================================================================

tau_', cov_clean, ' ~ T(
  dgamma(
    1,
    1
  ),
  0.00001,
  Inf
)

sigma_', cov_clean, ' <-
  1 / sqrt(tau_', cov_clean, ')
'
   )
    }

    # ------------------------------------------------------------------------
    # Biomass categorical random effect
    # ------------------------------------------------------------------------

    if (include_biomass) {

      random_effect_code <- c(
        random_effect_code,

  if(share_hyperparameters){
        paste0(
          '
# =============================================================================
# Biomass random effect: ', cov, '
# =============================================================================

for (j in 1:n_', cov_clean, ') {

  ', cov_clean, '_effect_biomass[j] ~ dnorm(
    0,
    tau_', cov_clean, '
  )

} '
   )
    
  } else {

paste0(
          '
# =============================================================================
# Biomass random effect: ', cov, '
# =============================================================================

for (j in 1:n_', cov_clean, ') {

  ', cov_clean, '_effect_biomass[j] ~ dnorm(
    0,
    tau_', cov_clean, '_biomass
  )

}

tau_', cov_clean, '_biomass ~ T(
  dgamma(
    1,
    1
  ),
  0.00001,
  Inf
)

sigma_', cov_clean, '_biomass <-
  1 / sqrt(tau_', cov_clean, '_biomass)
'
   )
  }
      )
  

biomass_random_terms <- c(
  biomass_random_terms,

  paste0(
    cov_clean,
    "_effect_biomass[",
    cov_clean,
    "_id[i]]"
  )
)

  }


  # ------------------------------------------------------------------------
  # Diversity categorical random effect
  # ------------------------------------------------------------------------

  if (include_diversity) {

    random_effect_code <- c(
      random_effect_code,

    if(share_hyperparameters) {
      paste0(
        '
# =============================================================================
# Diversity random effect: ', cov, '
# =============================================================================

for (j in 1:n_', cov_clean, ') {

  ', cov_clean, '_effect_diversity[j] ~ dnorm(
    0,
    tau_', cov_clean, '
  )
}

'
      )
    } else {
      paste0(
        '
# =============================================================================
# Diversity random effect: ', cov, '
# =============================================================================

for (j in 1:n_', cov_clean, ') {

  ', cov_clean, '_effect_diversity[j] ~ dnorm(
    0,
    tau_', cov_clean, '_diversity
  )

}

tau_', cov_clean, '_diversity ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_', cov_clean, '_diversity <-
  1 / sqrt(tau_', cov_clean, '_diversity)
'
      )
    }

    )


diversity_random_terms <- c(
  diversity_random_terms,

  paste0(
    cov_clean,
    "_effect_diversity[",
    cov_clean,
    "_id[i]]"
  )
)

  }

  }

}

# ============================================================================
# 15. Year random effect
# ============================================================================

if (random_year) {

  if(share_hyperparameters){
   
      '
# =============================================================================
# Biomass year random effect hyperameters
# =============================================================================

tau_year ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_year <-
  1 / sqrt(tau_year)
'
    
  }

  if (include_biomass) {

    random_effect_code <- c(
      random_effect_code,
      
    if(share_hyperparameters){
      '
# =============================================================================
# Biomass year random effect
# =============================================================================

for (j in 1:n_year) {

  year_effect_biomass[j] ~ dnorm(
    0,
    tau_year
  )

}
'
      
    } else {
      '
# =============================================================================
# Biomass year random effect
# =============================================================================

for (j in 1:n_year) {

  year_effect_biomass[j] ~ dnorm(
    0,
    tau_year_biomass
  )

}

tau_year_biomass ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_year_biomass <-
  1 / sqrt(tau_year_biomass)
'
    } 


    )


biomass_random_terms <- c(
  biomass_random_terms,

  "year_effect_biomass[year_id[i]]"
)

  }


if (include_diversity) {

  random_effect_code <- c(
    random_effect_code,

if(share_hyperparameters){

      '
# =============================================================================
# Diversity year random effect
# =============================================================================

for (j in 1:n_year) {

  year_effect_diversity[j] ~ dnorm(
    0,
    tau_year
  )

}
'
  
} else {

      '
# =============================================================================
# Diversity year random effect
# =============================================================================

for (j in 1:n_year) {

  year_effect_diversity[j] ~ dnorm(
    0,
    tau_year_diversity
  )

}

tau_year_diversity ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_year_diversity <-
  1 / sqrt(tau_year_diversity)
'
  
}
  )


diversity_random_terms <- c(
  diversity_random_terms,

  "year_effect_diversity[year_id[i]]"
)

}

}


# ============================================================================
# 16. Locality random effect
# ============================================================================

if (random_locality) {

  if(share_hyperparameters){
       '
# =============================================================================
# Locality hyperparameters
# =============================================================================
tau_locality ~ T(
  dgamma(
    1,
    1
  ),
  0.00001,
  Inf
)

sigma_locality <-
  1 / sqrt(tau_locality)
'   
  }

  if (include_biomass) {

    random_effect_code <- c(
      random_effect_code,

if(share_hyperparameters){

        '
# =============================================================================
# Biomass locality random effect
# =============================================================================

for (j in 1:n_locality) {

  locality_effect_biomass[j] ~ dnorm(
    0,
    tau_locality
  )

}

'
  
} else {

        '
# =============================================================================
# Biomass locality random effect
# =============================================================================

for (j in 1:n_locality) {

  locality_effect_biomass[j] ~ dnorm(
    0,
    tau_locality_biomass
  )

}

tau_locality_biomass ~ T(
  dgamma(
    1,
    1
  ),
  0.00001,
  Inf
)

sigma_locality_biomass <-
  1 / sqrt(tau_locality_biomass)
'

}
    )


biomass_random_terms <- c(
  biomass_random_terms,

  "locality_effect_biomass[locality_id[i]]"
)

  }


if (include_diversity) {

  random_effect_code <- c(
    random_effect_code,

    if(share_hyperparameters){

          '
# =============================================================================
# Diversity locality random effect
# =============================================================================

for (j in 1:n_locality) {

  locality_effect_diversity[j] ~ dnorm(
    0,
    tau_locality
  )
}
'

    } else {

          '
# =============================================================================
# Diversity locality random effect
# =============================================================================

for (j in 1:n_locality) {

  locality_effect_diversity[j] ~ dnorm(
    0,
    tau_locality_diversity
  )

}

tau_locality_diversity ~ T(
  dgamma(
    1,
    1
  ),
  0.00001,
  Inf
)

sigma_locality_diversity <-
  1 / sqrt(tau_locality_diversity)
'
    }


  )

diversity_random_terms <- c(
  diversity_random_terms,

  "locality_effect_diversity[locality_id[i]]"
)

}

}


# ============================================================================
# 17. Replicate random effect
# ============================================================================

if (random_replicate) {

  if(share_hyperparameters){
          '
# =============================================================================
# Replicate hyperameters
# =============================================================================

tau_replicate ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_replicate <-
  1 / sqrt(tau_replicate)
'
  }

  if (include_biomass) {

    random_effect_code <- c(
      random_effect_code,

if(share_hyperparameters){
      '
# =============================================================================
# Biomass replicate random effect
# =============================================================================

for (j in 1:n_replicate) {

  replicate_effect_biomass[j] ~ dnorm(
    0,
    tau_replicate
  )
}
'
} else {
  
      '
# =============================================================================
# Biomass replicate random effect
# =============================================================================

for (j in 1:n_replicate) {

  replicate_effect_biomass[j] ~ dnorm(
    0,
    tau_replicate_biomass
  )

}

tau_replicate_biomass ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_replicate_biomass <-
  1 / sqrt(tau_replicate_biomass)
'
  
}

    )


biomass_random_terms <- c(
  biomass_random_terms,

  "replicate_effect_biomass[replicate_id[i]]"
)

  }


if (include_diversity) {

  random_effect_code <- c(
    random_effect_code,

    if(share_hyperparameters){

        '
# =============================================================================
# Diversity replicate random effect
# =============================================================================

for (j in 1:n_replicate) {

  replicate_effect_diversity[j] ~ dnorm(
    0,
    tau_replicate
  )
}
'
      
    } else {

        '
# =============================================================================
# Diversity replicate random effect
# =============================================================================

for (j in 1:n_replicate) {

  replicate_effect_diversity[j] ~ dnorm(
    0,
    tau_replicate_diversity
  )

}

tau_replicate_diversity ~ T(
  dgamma(
    1,
    1
  ),
  0.00001,
  Inf
)

sigma_replicate_diversity <-
  1 / sqrt(tau_replicate_diversity)
'
      
    }
  
  )


diversity_random_terms <- c(
  diversity_random_terms,

  "replicate_effect_diversity[replicate_id[i]]"
)

}

}
  
# ===========================================================================0
# 18. Temporal effect
# ============================================================================

temporal_code <- character(0)

biomass_temporal_term <- character(0)

diversity_temporal_term <- character(0)



# ============================================================================
# 19. Biomass linear predictor
# ============================================================================

biomass_mean_terms <- c(
  "alpha_biomass",
  fixed_effect_biomass,
  biomass_random_terms,
  biomass_spatial_terms
)

biomass_mean_terms <- biomass_mean_terms[
  !is.na(biomass_mean_terms) &
    nzchar(biomass_mean_terms)
]

biomass_mean_expression <- paste(
  biomass_mean_terms,
  collapse = " + "
)


# ============================================================================
# 20. Diversity linear predictor
# ============================================================================

diversity_mean_terms <- c(
  "alpha_diversity",
  fixed_effect_diversity,
  diversity_random_terms,
  diversity_spatial_terms
)

diversity_mean_terms <- diversity_mean_terms[
  !is.na(diversity_mean_terms) &
    nzchar(diversity_mean_terms)
]

diversity_mean_expression <- paste(
  diversity_mean_terms,
  collapse = " + "
)


# ============================================================================
# 21. Biomass observation model
# ============================================================================

biomass_code <- ""

if (include_biomass) {

  if (biomass_distribution == "normal") {

    biomass_code <- paste0(

      '
# =============================================================================
# Biomass observation model
# Distribution: Normal
# =============================================================================

for (i in 1:N) {

  biomass[i] ~ dnorm(
    mu_biomass[i],
    tau_biomass
  )

  mu_biomass[i] <- ',
      biomass_mean_expression,
      '

}

# =============================================================================
# Biomass intercept
# =============================================================================

alpha_biomass ~ dnorm(
  0,
  0.001
)

# =============================================================================
# Biomass residual precision
# =============================================================================

tau_biomass ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_biomass <-
  1 / sqrt(tau_biomass)

'
    )

  }


if (biomass_distribution == "lognormal") {

  biomass_code <- paste0(

    '
# =============================================================================
# Biomass observation model
# Distribution: Log-normal
# =============================================================================

for (i in 1:N) {

  biomass[i] ~ dlnorm(
    mu_log_biomass[i],
    tau_biomass
  )

  mu_log_biomass[i] <- ',
    biomass_mean_expression,
    '

}

# =============================================================================
# Biomass intercept
# =============================================================================

alpha_biomass ~ dnorm(
  0,
  0.001
)

# =============================================================================
# Biomass residual precision
# =============================================================================

tau_biomass ~ T(
  dgamma(
    1,
    1
  ),
  0.0001,
  Inf
)

sigma_biomass <-
  1 / sqrt(tau_biomass)

'
  )

}

}


# ============================================================================
# 22. Diversity observation model
# ============================================================================

diversity_code <- ""

if (include_diversity) {

  if (diversity_distribution == "poisson") {

    diversity_code <- paste0(

      '
# =============================================================================
# Diversity observation model
# Distribution: Poisson
# =============================================================================

for (i in 1:N) {

  diversity[i] ~ dpois(
    lambda_diversity[i]
  )

  log(lambda_diversity[i]) <- ',
      diversity_mean_expression,
      '

}

# =============================================================================
# Diversity intercept
# =============================================================================

alpha_diversity ~ dnorm(
  0,
  0.001
)

'
    )

  }


if (diversity_distribution == "negative_binomial") {

  diversity_code <- paste0(

    '
# =============================================================================
# Diversity observation model
# Distribution: Negative binomial
# =============================================================================

for (i in 1:N) {

  diversity[i] ~ dnbinom(
    size_diversity,
    prob_diversity[i]
  )

  prob_diversity[i] <-
    size_diversity /
    (size_diversity + mu_diversity[i])

  log(mu_diversity[i]) <- ',
    diversity_mean_expression,
    '

}

# =============================================================================
# Diversity intercept
# =============================================================================

alpha_diversity ~ dnorm(
  0,
  0.001
)

# =============================================================================
# Negative-binomial dispersion
# =============================================================================

size_diversity ~ dgamma(
  0.001,
  0.001
)

'
  )

}

}


# ============================================================================
# 23. Shared fixed-effect priors
# ============================================================================

beta_code <- character(0)

if (n_continuous > 0) {

  if(share_fixed_effects){

      beta_code <- c(
    beta_code,

    '
# =============================================================================
# Diversity continuous covariate effects
# =============================================================================

for (j in 1:p) {

  beta[j] ~ dnorm(
    0,
    0.001
  )

}
'
  )

  } else {

      if (include_biomass) {

    beta_code <- c(
      beta_code,

      '
# =============================================================================
# Biomass continuous covariate effects
# =============================================================================

for (j in 1:p) {

  beta_biomass[j] ~ dnorm(
    0,
    0.001
  )

}
'
    )

  }


if (include_diversity) {

  beta_code <- c(
    beta_code,

    '
# =============================================================================
# Diversity continuous covariate effects
# =============================================================================

for (j in 1:p) {

  beta_diversity[j] ~ dnorm(
    0,
    0.001
  )

}
'
  )

}
    
  }

}


# 24. Spatial hyperparameter constants
# ============================================================================

spatial_hyperparameter_code <- character(0)

if (spatial_effect == "autocorrelation") {

  spatial_hyperparameter_code <- paste0(

    '
# =============================================================================
# Spatial prior hyperparameters
# =============================================================================

spatial_range_lower <- ',
    spatial_range_prior[1],
    '

spatial_range_upper <- ',
    spatial_range_prior[2],
    '

spatial_sigma_precision <- ',
    1 / spatial_sigma_prior^2,
    '

'
  )

}


# ============================================================================
# 25. Construct complete NIMBLE model
# ============================================================================

nimble_model <- paste0(

  '
# =============================================================================
# NIMBLE model: insect biomass and diversity
# =============================================================================

# Response:
#   ', response, '

# Biomass distribution:
#   ', biomass_distribution, '

# Diversity distribution:
#   ', diversity_distribution, '

# Spatial effect:
#   ', spatial_effect, '

# Shared spatial process:
#   ', spatial_shared, '

# Effect structure:
#   ', effect_structure, '

# Temporal effect:
#   ', temporal_effect, '

# Shared temporal process:
#   ', temporal_shared, '

# Continuous covariates:
#   ',
  ifelse(
    length(continuous_covariates) > 0,
    paste(
      continuous_covariates,
      collapse = ", "
    ),
    "None"
  ),
  '

# Categorical covariates:
#   ',
  ifelse(
    length(categorical_covariates) > 0,
    paste(
      categorical_covariates,
      collapse = ", "
    ),
    "None"
  ),
  '

# =============================================================================

library(nimble)

biomass_diversity_model <- nimbleCode({

',
  spatial_hyperparameter_code,
  "\n",

  biomass_code,
  diversity_code,
  "\n",

  paste(
    beta_code,
    collapse = "\n"
  ),
  "\n",

  paste(
    random_effect_code,
    collapse = "\n"
  ),
  "\n",

  paste(
    spatial_code,
    collapse = "\n"
  ),
  "\n",

  '
})

'
)


# ============================================================================
# 26. Create output directory
# ============================================================================

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================================
# 27. Model output file
# ============================================================================

nimble_file <- file.path(
  output_dir,
  paste0(
    model_name,
    ".R"
  )
)


# ============================================================================
# 28. Write NIMBLE model to file
# ============================================================================

writeLines(
  nimble_model,
  con = nimble_file
)


# ============================================================================
# 29. Return model information
# ============================================================================

invisible(
  list(

    # ------------------------------------------------------------------------
    # Generated model
    # ------------------------------------------------------------------------

    model = nimble_model,

    file = nimble_file,


    # ------------------------------------------------------------------------
    # Response specification
    # ------------------------------------------------------------------------

    response = response,


    # ------------------------------------------------------------------------
    # Distribution specification
    # ------------------------------------------------------------------------

    biomass_distribution = biomass_distribution,

    diversity_distribution = diversity_distribution,


    # ------------------------------------------------------------------------
    # Covariates
    # ------------------------------------------------------------------------

    covariates = modelCovs,

    continuous_covariates = continuous_covariates,

    categorical_covariates = categorical_covariates,

    covariate_types = covariate_types,


    # ------------------------------------------------------------------------
    # Random effects
    # ------------------------------------------------------------------------

    random_year = random_year,

    random_locality = random_locality,

    random_replicate = random_replicate,


    # ------------------------------------------------------------------------
    # Temporal and spatio-temporal structure
    # ------------------------------------------------------------------------

    effect_structure = effect_structure,

    temporal_effect = temporal_effect,

    temporal_shared = temporal_shared,

    temporal_basis = temporal_basis,

    temporal_penalty = temporal_penalty,

    temporal_basis_penalized = temporal_basis_penalized,

    temporal_penalty_eigenvalues =
      temporal_penalty_eigenvalues,

    temporal_penalty_eigenvectors =
      temporal_penalty_eigenvectors,

    temporal_penalty_tol =
      temporal_penalty_tol,

    n_temporal_basis =
      n_temporal_basis,

    # These objects should be supplied as NIMBLE constants
    temporal_constants = list(

      n_temporal_basis =
        n_temporal_basis,

      temporal_basis_penalized =
        temporal_basis_penalized,

      temporal_penalty_eigenvalues =
        temporal_penalty_eigenvalues

    ),


    # ------------------------------------------------------------------------
    # Spatial effect
    # ------------------------------------------------------------------------

    spatial_effect = spatial_effect,

    # The current spatiotemporal option combines additive spatial and temporal
    # effects. A location-specific temporal interaction requires an additional
    # space-by-time basis or dynamic spatial process and is not generated here.

    spatial_shared = spatial_shared,

    spatial_range_prior = spatial_range_prior,

    spatial_sigma_prior = spatial_sigma_prior

  )
)

}
