#' Create a Temporal Insect Biomass and Diversity Model
#'
#' Generates an R/NIMBLE model specification for insect biomass, diversity,
#' or both responses with a temporal component. The generated model code is
#' written to an R file, and the function invisibly returns the code, output
#' path, and model configuration.
#'
#' The temporal component can be a quadratic linear-time effect or a penalized
#' cubic cyclic spline. By default, the spline process is response-specific.
#' For cubic cyclic splines, `temporal_shared = TRUE` allows the responses to
#' use a shared latent temporal curve, with a response-specific loading when
#' both responses are included.
#'
#' @param modelCovs Character vector naming the covariates to include as
#'   predictors. Corresponding design-matrix inputs and any required
#'   covariate-level indices must be supplied when constructing the model.
#' @param response Response to model: `"biomass"`, `"diversity"`, or `"both"`.
#'   Defaults to `"biomass"`.
#' @param biomass_distribution Biomass distribution, either `"normal"` or
#'   `"lognormal"`. Defaults to `"normal"`.
#' @param diversity_distribution Diversity distribution, either `"poisson"` or
#'   `"negative_binomial"`. Defaults to `"poisson"`.
#' @param covariate_types Optional named character vector assigning each
#'   covariate to `"continuous"` or `"categorical"`. If omitted, all covariates
#'   are treated as continuous.
#' @param random_year Logical; include response-specific year random effects
#'   when `TRUE`. Defaults to `FALSE`.
#' @param random_locality Logical; include response-specific locality random
#'   effects when `TRUE`. Defaults to `FALSE`.
#' @param random_replicate Logical; include response-specific replicate random
#'   effects when `TRUE`. Defaults to `FALSE`.
#' @param temporal_effect Temporal effect: `"none"`, `"linear"`, or
#'   `"cubic_cyclic"`. The cyclic spline requires `temporal_basis` and
#'   `temporal_penalty`. Defaults to `"cubic_cyclic"`.
#' @param temporal_basis Numeric matrix of temporal basis evaluations, required
#'   for `"cubic_cyclic"`, with rows corresponding to observations.
#' @param temporal_penalty Numeric square penalty matrix with dimensions
#'   matching the number of columns in `temporal_basis`; required for
#'   `"cubic_cyclic"`.
#' @param temporal_penalty_tol Positive finite tolerance used relative to the
#'   largest penalty eigenvalue to retain the positive-eigenvalue subspace.
#'   Defaults to `1e-8`.
#' @param temporal_shared Logical; for `"cubic_cyclic"`, use a shared latent
#'   temporal process if `TRUE`. If both responses are included, the biomass
#'   loading is fixed to one and the diversity loading is estimated. With
#'   `FALSE`, response-specific spline coefficients and smoothing precisions
#'   are used. Shared temporal effects are not supported for `"linear"`.
#' @param spatial_effect Spatial-effect option accepted for interface
#'   consistency. Spatial effects are disabled in this temporal-only builder.
#' @param spatial_shared Spatial-process sharing option for related builders;
#'   it does not add a spatial process to this temporal-only model.
#' @param spatial_range_prior Lower and upper bounds for the spatial range
#'   prior in spatial model builders. Not used by this temporal-only model.
#' @param spatial_sigma_prior Positive scale for the spatial standard-deviation
#'   prior in spatial model builders. Not used by this temporal-only model.
#' @param output_dir Directory in which to write the generated model file.
#'   Defaults to `temp_dir`.
#' @param model_name Base name of the generated R file, without `.R`.
#' @details
#' The model has response-specific intercepts and fixed-effect coefficients
#' for continuous covariates. Categorical covariates and optional year,
#' locality, and replicate effects are response-specific.
#'
#' For `"linear"`, the temporal contribution includes first- and second-order
#' terms in `temporal[i]`. For `"cubic_cyclic"`, the basis is rotated into the
#' positive-eigenvalue basis of the penalty matrix; coefficient prior
#' precisions are scaled by retained penalty eigenvalues and a smoothing
#' precision parameter.
#'
#' This is a temporal-only model builder. It generates NIMBLE code but does
#' not construct, compile, or fit the model. Compatible data, constants, and
#' initial values must be supplied separately.
#' @return Invisibly returns a list containing generated model code and file
#'   path, response and distribution settings, covariate specifications,
#'   random-effect settings, temporal basis/penalty information, and model
#'   configuration metadata.
#' @seealso
#' \code{\link{create_insect_spatiotemporal_shared_model}},
#' \code{\link{create_insect_spatial_separate_model}},
#' \code{\link{create_insect_nimble_model}}
#' @examples
#' \dontrun{
#' model_info <- create_insect_temporal_separate_model(
#'   modelCovs = c("temperature", "vegetation"),
#'   response = "both",
#'   covariate_types = c(
#'     temperature = "continuous",
#'     vegetation = "categorical"
#'   ),
#'   temporal_effect = "cubic_cyclic",
#'   temporal_basis = temporal_basis,
#'   temporal_penalty = temporal_penalty,
#'   temporal_shared = FALSE,
#'   output_dir = "model_code"
#' )
#' cat(model_info$model)
#' model_info$file
#' }
#' @export


create_insect_temporal_separate_model <- function(

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

  temporal_effect = c("cubic_cyclic", "linear"),

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

  spatial_effect = "none",

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

  model_name = "create_temporal_separate_model"

) {

  # Fixed structure for this standalone model builder: Temporal-only model with separate response-specific temporal parameters.
  effect_structure <- "temporal"


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

  if (temporal_shared && temporal_effect != "cubic_cyclic") {
    stop(
      "'temporal_shared = TRUE' is currently supported only with ",
      "temporal_effect = 'cubic_cyclic'.",
      call. = FALSE
    )
  }

  # Validate the relative eigenvalue threshold used to identify
  # the positive-eigenvalue subspace of the spline penalty.
  if (
    !is.numeric(temporal_penalty_tol) ||
    length(temporal_penalty_tol) != 1L ||
    !is.finite(temporal_penalty_tol) ||
    temporal_penalty_tol <= 0
  ) {
    stop(
      "'temporal_penalty_tol' must be a single positive, finite number.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------------------------
  # Validate and prepare cubic cyclic spline
  # ----------------------------------------------------------------------------

  if (temporal_effect == "cubic_cyclic") {

    if (is.null(temporal_basis)) {

      stop(
        "'temporal_basis' must be supplied when ",
        "temporal_effect = 'cubic_cyclic'.",
        call. = FALSE
      )

    }

    if (is.null(temporal_penalty)) {

      stop(
        "'temporal_penalty' must be supplied when ",
        "temporal_effect = 'cubic_cyclic'.",
        call. = FALSE
      )

    }

    if (
      !is.matrix(temporal_basis) ||
      !is.numeric(temporal_basis)
    ) {

      stop(
        "'temporal_basis' must be a numeric matrix.",
        call. = FALSE
      )

    }

    if (
      !is.matrix(temporal_penalty) ||
      !is.numeric(temporal_penalty)
    ) {

      stop(
        "'temporal_penalty' must be a numeric matrix.",
        call. = FALSE
      )

    }

    n_basis_original <- ncol(temporal_basis)

    if (
      nrow(temporal_penalty) != n_basis_original ||
      ncol(temporal_penalty) != n_basis_original
    ) {

      stop(
        "'temporal_penalty' must have dimensions ",
        n_basis_original,
        " x ",
        n_basis_original,
        ".",
        call. = FALSE
      )

    }

    if (
      any(!is.finite(temporal_basis)) ||
      any(!is.finite(temporal_penalty))
    ) {

      stop(
        "'temporal_basis' and 'temporal_penalty' ",
        "must contain only finite values.",
        call. = FALSE
      )

    }


    # --------------------------------------------------------------------------
    # Force numerical symmetry
    # --------------------------------------------------------------------------

    temporal_penalty <- (
      temporal_penalty +
        t(temporal_penalty)
    ) / 2


    # --------------------------------------------------------------------------
    # Eigen-decomposition of penalty matrix
    #
    # S = V D V'
    #
    # The cyclic cubic spline penalty contains a null space.
    # We retain only the positive eigenvalues.
    # --------------------------------------------------------------------------

    temporal_eigen <- eigen(
      temporal_penalty,
      symmetric = TRUE
    )

    temporal_scale <- max(
      temporal_eigen$values
    )

    if (
      !is.finite(temporal_scale) ||
      temporal_scale <= 0
    ) {

      stop(
        "'temporal_penalty' must contain at least one ",
        "positive eigenvalue.",
        call. = FALSE
      )

    }


    temporal_keep <- (
      temporal_eigen$values >
        temporal_penalty_tol * temporal_scale
    )


    if (!any(temporal_keep)) {

      stop(
        "No positive eigenvalues remain in 'temporal_penalty'.",
        call. = FALSE
      )

    }


    # --------------------------------------------------------------------------
    # Penalized eigenvalues and eigenvectors
    # --------------------------------------------------------------------------

    temporal_penalty_eigenvalues <- (
      temporal_eigen$values[
        temporal_keep
      ]
    )

    temporal_penalty_eigenvectors <- (
      temporal_eigen$vectors[
        ,
        temporal_keep,
        drop = FALSE
      ]
    )


    # --------------------------------------------------------------------------
    # Rotate basis into penalized eigenbasis
    #
    # B_pen = B V_+
    #
    # This allows the coefficient prior to be represented as independent
    # normals with precision tau * eigenvalue.
    # --------------------------------------------------------------------------

    temporal_basis_penalized <- (
      temporal_basis %*%
        temporal_penalty_eigenvectors
    )

    n_temporal_basis <- ncol(
      temporal_basis_penalized
    )


  } else {

    temporal_basis_penalized <- NULL

    temporal_penalty_eigenvalues <- NULL

    temporal_penalty_eigenvectors <- NULL

    n_temporal_basis <- 0L

  }


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



  # ============================================================================
  # 10. Initialise random effects
  # ============================================================================

  random_effect_code <- character(0)

  biomass_random_terms <- character(0)

  diversity_random_terms <- character(0)

# ============================================================================
# 14. Categorical covariates as random effects
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

# ============================================================================
# 18. Temporal effect
# ============================================================================

temporal_code <- character(0)

biomass_temporal_term <- character(0)

diversity_temporal_term <- character(0)


# ----------------------------------------------------------------------------
# Linear temporal effect
# ----------------------------------------------------------------------------

if (temporal_effect == "linear") {

  if (include_biomass) {

    temporal_code <- c(
      temporal_code,

      '
# =============================================================================
# Biomass linear temporal effect
# =============================================================================

beta_time_biomass ~ dnorm(
  0,
  0.001
)

beta_time_biomass_squared ~ dnorm(
  0,
  0.001
)
'
    )


biomass_temporal_term <- c(
  biomass_temporal_term,

  "beta_time_biomass * temporal[i]",

  "beta_time_biomass_squared * pow(temporal[i], 2)"
)

  }


  if (include_diversity) {

    temporal_code <- c(
      temporal_code,

      '
# =============================================================================
# Diversity linear temporal effect
# =============================================================================

beta_time_diversity ~ dnorm(
  0,
  0.001
)

beta_time_diversity_squared ~ dnorm(
  0,
  0.001
)
'
    )


diversity_temporal_term <- c(
  diversity_temporal_term,

  "beta_time_diversity * temporal[i]",

  "beta_time_diversity_squared * pow(temporal[i], 2)"
)

  }

}


# ----------------------------------------------------------------------------
# Penalized cubic cyclic spline
# ----------------------------------------------------------------------------

if (temporal_effect == "cubic_cyclic") {

  if (temporal_shared) {

    # A single latent temporal curve is shared across the response models.
    # For both responses, fix the biomass loading to 1 to identify the scale;
    # estimate only the diversity loading. With one response, its loading is 1.
    temporal_code <- c(
      temporal_code,
      '
# =============================================================================
# Shared penalized cubic cyclic temporal process
# =============================================================================

for (k in 1:n_temporal_basis) {
  beta_time_shared[k] ~ dnorm(
    0,
    temporal_tau_shared * temporal_penalty_eigenvalues[k]
  )
}

temporal_tau_shared ~ dgamma(1, 1)
'
    )

    shared_temporal_expression <- paste0(
      "inprod(beta_time_shared[1:n_temporal_basis], ",
      "temporal_basis_penalized[i, 1:n_temporal_basis])"
    )

    if (include_biomass) {
      biomass_temporal_term <- c(
        biomass_temporal_term,
        shared_temporal_expression
      )
    }

    if (include_diversity) {
      if (include_biomass) {
        temporal_code <- c(
          temporal_code,
          '
# Biomass loading fixed to 1 for identifiability;
# estimate the diversity-specific loading.
lambda_time_diversity ~ dnorm(0, 1)
'
        )
        diversity_temporal_term <- c(
          diversity_temporal_term,
          paste0("lambda_time_diversity * (", shared_temporal_expression, ")")
        )
      } else {
        # Single-response model: the loading is fixed to 1.
        diversity_temporal_term <- c(
          diversity_temporal_term,
          shared_temporal_expression
        )
      }
    }

  } else {

    if (include_biomass) {
      temporal_code <- c(
        temporal_code,
        '
# =============================================================================
# Penalized cubic cyclic spline: biomass
# =============================================================================

for (k in 1:n_temporal_basis) {
  beta_time_biomass[k] ~ dnorm(
    0,
    temporal_tau_biomass * temporal_penalty_eigenvalues[k]
  )
}

temporal_tau_biomass ~ dgamma(1, 1)
'
      )

      biomass_temporal_term <- c(
        biomass_temporal_term,
        paste0(
          "inprod(beta_time_biomass[1:n_temporal_basis], ",
          "temporal_basis_penalized[i, 1:n_temporal_basis])"
        )
      )
    }

if (include_diversity) {
  temporal_code <- c(
    temporal_code,
    '
# =============================================================================
# Penalized cubic cyclic spline: diversity
# =============================================================================

for (k in 1:n_temporal_basis) {
  beta_time_diversity[k] ~ dnorm(
    0,
    temporal_tau_diversity * temporal_penalty_eigenvalues[k]
  )
}

temporal_tau_diversity ~ dgamma(1, 1)
'
  )

  diversity_temporal_term <- c(
    diversity_temporal_term,
    paste0(
      "inprod(beta_time_diversity[1:n_temporal_basis], ",
      "temporal_basis_penalized[i, 1:n_temporal_basis])"
    )
  )
}
  }
}


# ============================================================================
# 19. Biomass linear predictor
# ============================================================================

biomass_mean_terms <- c(
  "alpha_biomass",
  fixed_effect_biomass,
  biomass_random_terms,
  biomass_temporal_term
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
  diversity_temporal_term
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

  biomass[i] ~ T(dnorm(
    mu_biomass[i],
    tau_biomass
  ), 0, )

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
# 23. Fixed-effect priors
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

# ============================================================================
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
    temporal_code,
    collapse = "\n"
  ),

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
