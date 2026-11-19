#' Create a Baseline NIMBLE Model for Insect Biomass and Diversity
#'
#' Generate a standalone NIMBLE model for insect biomass, diversity, or both
#' responses without temporal or spatial processes. The model supports
#' alternative response distributions, continuous covariate effects,
#' categorical covariates represented as random effects, and optional
#' random intercepts for year, locality, and replicate.
#'
#' The function constructs the NIMBLE model code as text, writes it to an R
#' script, and invisibly returns the generated code together with metadata
#' describing the model specification. The resulting script can subsequently
#' be sourced to define the NIMBLE model.
#'
#' For biomass, the observation distribution can be normal or log-normal.
#' For diversity, the observation distribution can be Poisson or negative
#' binomial. When both responses are requested, separate response-specific
#' linear predictors and regression coefficients are generated.
#'
#' Continuous covariates enter the linear predictors as fixed effects.
#' Categorical covariates are represented by level-specific random effects
#' with normally distributed effects and estimated precision parameters.
#' Optional year, locality, and replicate effects are also represented by
#' response-specific random intercepts.
#'
#' This function is intended to provide a baseline candidate model against
#' which models containing temporal, spatial, or spatiotemporal processes
#' can be compared. The baseline model does not require temporal basis
#' matrices, spatial coordinates, or spatial autocorrelation structures.
#'
#' @param modelCovs Character vector giving the names of covariates to include
#'   in the model. Covariates classified as continuous are included as fixed
#'   effects; covariates classified as categorical are represented as
#'   random effects. Use `character(0)` when no covariates are required.
#'
#' @param response Character string specifying which response to model.
#'   Must be one of `"biomass"`, `"diversity"`, or `"both"`.
#'   The default is `"biomass"`.
#'
#' @param biomass_distribution Character string specifying the biomass
#'   observation distribution. Must be `"normal"` or `"lognormal"`.
#'   The default is `"normal"`. This argument is relevant when biomass
#'   is included in the model.
#'
#' @param diversity_distribution Character string specifying the diversity
#'   observation distribution. Must be `"poisson"` or
#'   `"negative_binomial"`. The default is `"poisson"`. This argument is
#'   relevant when diversity is included in the model.
#'
#' @param covariate_types Optional named character vector identifying the
#'   type of each covariate in `modelCovs`. Names must correspond to
#'   covariate names, and values must be `"continuous"` or `"categorical"`.
#'   If `NULL`, all specified covariates are treated as continuous.
#'   Every covariate in `modelCovs` must have a corresponding type when
#'   this argument is supplied.
#' 
#' @param share_fixed_effects Logical. If `TRUE`, all the fixed effects will
#' be shared by both the response variables. Defaults to `FALSE`.
#' 
#' @param share_hyperpameters Logical. If `TRUE`, all random effect hyperparameters
#' will be shared by both response variables. Defauls to `FALSE`.
#'
#' @param random_year Logical. If `TRUE`, include a year-specific random
#'   intercept for each response included in the model. Defaults to
#'   `FALSE`.
#'
#' @param random_locality Logical. If `TRUE`, include a locality-specific
#'   random intercept for each response included in the model. Defaults
#'   to `FALSE`.
#'
#' @param random_replicate Logical. If `TRUE`, include a replicate-specific
#'   random intercept for each response included in the model. Defaults
#'   to `FALSE`.
#'
#' @param temporal_effect Character string specifying the temporal effect.
#'   The baseline model uses `"none"`. The argument is retained for
#'   compatibility with the model-building interface, but this function
#'   is intended to generate a model without temporal processes.
#'
#' @param temporal_basis Optional matrix containing a temporal spline basis.
#'   Not required for the baseline model.
#'
#' @param temporal_penalty Optional matrix containing the penalty matrix
#'   associated with a temporal spline basis. Not required for the baseline
#'   model.
#'
#' @param temporal_penalty_tol Positive finite numeric value specifying
#'   the relative eigenvalue threshold used to identify the penalized
#'   subspace of a temporal spline penalty. Defaults to `1e-8`.
#'   Not required for the baseline model.
#'
#' @param temporal_shared Logical. Intended to indicate whether a temporal
#'   process is shared between the two responses. Defaults to `FALSE`.
#'   The baseline model does not include a temporal process.
#'
#' @param spatial_effect Character string specifying the spatial effect.
#'   The baseline model uses `"none"`. The argument is retained for
#'   compatibility with the model-building interface, but this function
#'   is intended to generate a model without spatial processes.
#'
#' @param spatial_shared Logical. Indicates whether a spatial
#'   autocorrelation process is shared between responses in model
#'   structures that support such a process. Defaults to `TRUE`.
#'   This option is not relevant to the baseline model.
#'
#' @param spatial_range_prior Numeric vector of length two specifying the
#'   lower and upper bounds used for the spatial range prior in model
#'   structures that support spatial autocorrelation. Both values must be
#'   positive, and the upper bound must exceed the lower bound. Defaults
#'   to `c(1000, 500000)`. This argument is not relevant to the baseline
#'   model.
#'
#' @param spatial_sigma_prior Positive finite numeric value specifying
#'   the spatial standard-deviation prior parameter in model structures
#'   that support spatial autocorrelation. Defaults to `1`. This argument
#'   is not relevant to the baseline model.
#'
#' @param output_dir Character string specifying the directory in which
#'   to write the generated model script. Defaults to `temp_dir`.
#'   The directory is created recursively if it does not already exist.
#'
#' @param model_name Character string specifying the base name of the
#'   generated model script, without the `.R` extension. Defaults to
#'   `"create_insect_nimble_baseline_model"`.
#'
#' @return Invisibly returns a named list containing:
#' \describe{
#'   \item{model}{Character string containing the generated NIMBLE model code.}
#'   \item{file}{Path to the generated R script.}
#'   \item{response}{Selected response specification.}
#'   \item{biomass_distribution}{Selected biomass distribution.}
#'   \item{diversity_distribution}{Selected diversity distribution.}
#'   \item{covariates}{Names of all specified covariates.}
#'   \item{continuous_covariates}{Names of covariates included as fixed effects.}
#'   \item{categorical_covariates}{Names of covariates represented as random effects.}
#'   \item{covariate_types}{Named vector specifying covariate types.}
#'   \item{random_year}{Whether year-level random effects were requested.}
#'   \item{random_locality}{Whether locality-level random effects were requested.}
#'   \item{random_replicate}{Whether replicate-level random effects were requested.}
#'   \item{effect_structure}{The effect-structure metadata returned by the function.}
#'   \item{temporal_effect}{The temporal-effect setting.}
#'   \item{temporal_shared}{Whether a shared temporal process was requested.}
#'   \item{temporal_basis}{The supplied temporal basis, if any.}
#'   \item{temporal_penalty}{The supplied temporal penalty matrix, if any.}
#'   \item{temporal_basis_penalized}{The penalized temporal basis, if generated.}
#'   \item{temporal_penalty_eigenvalues}{Eigenvalues associated with the temporal penalty, if generated.}
#'   \item{temporal_penalty_eigenvectors}{Eigenvectors associated with the temporal penalty, if generated.}
#'   \item{temporal_penalty_tol}{The temporal penalty eigenvalue tolerance.}
#'   \item{n_temporal_basis}{The number of temporal basis functions retained.}
#'   \item{temporal_constants}{A list of temporal quantities intended for use as NIMBLE constants.}
#'   \item{spatial_effect}{The spatial-effect setting.}
#'   \item{spatial_shared}{Whether a shared spatial process was requested.}
#'   \item{spatial_range_prior}{The spatial range prior specification.}
#'   \item{spatial_sigma_prior}{The spatial standard-deviation prior parameter.}
#' }
#'
#' @details
#' The generated biomass model uses an identity link for a normal response
#' and a log link for a log-normal response parameterized by its log-scale
#' location. The diversity model uses a log link for the Poisson mean or
#' negative-binomial mean.
#'
#' Intercepts and continuous covariate coefficients have normal priors.
#' Random effects are normally distributed with precision parameters
#' assigned gamma-based priors, subject to the truncation specified in the
#' generated model code. The negative-binomial model includes a
#' gamma-prior dispersion parameter.
#'
#' When `response = "both"`, the biomass and diversity models have separate
#' intercepts, continuous covariate coefficients, and response-specific
#' random effects. This function does not, by itself, impose a shared
#' latent process between the responses.
#'
#' The generated model assumes that the data list supplied when building
#' the NIMBLE model contains the response variables and any covariate,
#' grouping-index, and dimension objects required by the selected model
#' specification. Categorical covariates and optional grouping effects
#' require the corresponding level counts and observation-level indices.
#'
#' @examples
#' \dontrun{
#' # Baseline biomass model with continuous covariates
#' baseline_biomass <- create_insect_baseline_model(
#'   modelCovs = c("temperature", "precipitation"),
#'   response = "biomass",
#'   biomass_distribution = "normal",
#'   covariate_types = c(
#'     temperature = "continuous",
#'     precipitation = "continuous"
#'   ),
#'   output_dir = "models",
#'   model_name = "baseline_biomass"
#' )
#'
#' # Baseline model for both responses with grouping-level random effects
#' baseline_both <- create_insect_baseline_model(
#'   modelCovs = c("temperature", "habitat"),
#'   response = "both",
#'   biomass_distribution = "normal",
#'   diversity_distribution = "negative_binomial",
#'   covariate_types = c(
#'     temperature = "continuous",
#'     habitat = "categorical"
#'   ),
#'   random_year = TRUE,
#'   random_locality = TRUE,
#'   share_fixed_effects = FALSE,
#'   output_dir = "models",
#'   model_name = "baseline_biomass_diversity"
#' )
#' }
#'
#' @seealso
#' \code{\link[nimble]{nimbleCode}},
#' \code{\link[nimble]{configureMCMC}},
#' \code{\link[nimble]{buildMCMC}}
#'
#' @export

create_insect_baseline_model <- function(

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

  model_name = "create_insect_nimble_baseline_model"

) {

  # Fixed structure for this standalone model builder: Baseline model without temporal or spatial effects.
  effect_structure <- "none"


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
  # 11. Initialise spatial effects
  # ============================================================================

  spatial_code <- character(0)

  biomass_spatial_terms <- character(0)

  diversity_spatial_terms <- character(0)


  # ============================================================================
  # 12. Spatial coordinate effects
  # ============================================================================

  # ============================================================================
  # 13. Spatial autocorrelation
  # ============================================================================


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

# ----------------------------------------------------------------------------
# Penalized cubic cyclic spline
# ----------------------------------------------------------------------------

# ============================================================================
# 19. Biomass linear predictor
# ============================================================================

biomass_mean_terms <- c(
  "alpha_biomass",
  fixed_effect_biomass,
  biomass_random_terms,
  biomass_spatial_terms,
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
  diversity_spatial_terms,
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

