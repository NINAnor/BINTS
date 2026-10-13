#' Create a NIMBLE Model for Insect Biomass and Diversity
#'
#' Generate NIMBLE model code for Bayesian analyses of insect biomass,
#' diversity, or both responses. The function supports alternative response
#' distributions, continuous fixed effects, categorical random effects,
#' optional grouping-level random effects, temporal effects, and spatial
#' effects.
#'
#' Continuous covariates are included as fixed effects, whereas categorical
#' covariates are represented as random effects. Optional random intercepts
#' can be specified for year, locality, and replicate. Temporal structure can
#' be modelled using no temporal effect, a linear temporal effect, or a
#' penalized cubic cyclic spline. Spatial structure can be represented by
#' coordinate-based trends or a latent Gaussian spatial autocorrelation
#' process.
#'
#' For cyclic cubic splines, the function accepts a basis matrix and its
#' associated penalty matrix, typically obtained from `mgcv`. The penalty
#' matrix is symmetrized and decomposed into eigenvalues and eigenvectors.
#' Eigenvectors associated with positive eigenvalues are retained to construct
#' a penalized basis with independent Gaussian coefficient priors.
#'
#' The function constructs the model using `nimbleCode()` syntax, writes the
#' generated model code to an R script, and invisibly returns the model code,
#' file path, model configuration, covariate information, and relevant
#' temporal and spatial metadata.
#'
#' @param modelCovs A character vector specifying the names of covariates
#'   included in the model. Covariates must be represented in the model data
#'   using the corresponding continuous design-matrix columns or categorical
#'   index variables.
#'
#' @param response A character string specifying which response to model.
#'   Supported values are `"biomass"`, `"diversity"`, and `"both"`.
#'
#' @param biomass_distribution A character string specifying the observation
#'   distribution for biomass. Supported values are `"normal"` and
#'   `"lognormal"`. The default is `"normal"`.
#'
#' @param diversity_distribution A character string specifying the observation
#'   distribution for diversity. Supported values are `"poisson"` and
#'   `"negative_binomial"`. The default is `"poisson"`.
#'
#' @param covariate_types An optional named character vector identifying each
#'   covariate as `"continuous"` or `"categorical"`. Names must correspond to
#'   elements of `modelCovs`. If `NULL`, all covariates are treated as
#'   continuous.
#'
#' @param random_year Logical. If `TRUE`, include a year-level random
#'   intercept for each included response. Defaults to `FALSE`.
#'
#' @param random_locality Logical. If `TRUE`, include a locality-level random
#'   intercept for each included response. Defaults to `FALSE`.
#'
#' @param random_replicate Logical. If `TRUE`, include a replicate-level random
#'   intercept for each included response. Defaults to `FALSE`.
#'
#' @param temporal_effect A character string specifying the temporal
#'   structure. Supported values are `"none"`, `"linear"`, and
#'   `"cubic_cyclic"`. The default is `"none"`.
#'
#' @param temporal_basis A numeric matrix containing the temporal spline basis
#'   evaluated at the observations. Required when
#'   `temporal_effect = "cubic_cyclic"`. The number of columns determines the
#'   dimension of the original spline coefficient space.
#'
#' @param temporal_penalty A numeric square matrix representing the penalty
#'   associated with `temporal_basis`. Required when
#'   `temporal_effect = "cubic_cyclic"`. Its dimensions must equal the number
#'   of columns in `temporal_basis`.
#'
#' @param temporal_penalty_tol A positive numeric tolerance used to determine
#'   which eigenvalues of the temporal penalty matrix are retained. Eigenvalues
#'   are retained when they exceed this tolerance multiplied by the maximum
#'   eigenvalue. Defaults to `1e-8`.
#'
#' @param spatial_effect A character string specifying the spatial structure.
#'   Supported values are `"none"`, `"coordinates"`, and `"autocorrelation"`.
#'   The default is `"none"`.
#'
#' @param spatial_shared Logical. For `spatial_effect = "autocorrelation"`,
#'   specifies whether biomass and diversity share a latent spatial process.
#'   If `TRUE`, the process is shared and each response has its own loading.
#'   If `FALSE`, separate latent spatial processes are generated for the two
#'   responses. Defaults to `TRUE`.
#'
#' @param spatial_range_prior A numeric vector of length two specifying the
#'   lower and upper bounds of the uniform prior for the spatial range
#'   parameter. Both values must be positive, and the upper bound must exceed
#'   the lower bound. Defaults to `c(1000, 500000)`.
#'
#' @param spatial_sigma_prior A positive finite numeric scalar specifying the
#'   scale used to parameterize the half-normal prior for the spatial standard
#'   deviation. Defaults to `1`.
#'
#' @param output_dir A character string specifying the directory in which the
#'   generated model script is saved. The directory is created if it does not
#'   exist. Defaults to `temp_dir`.
#'
#' @param model_name A character string specifying the base name of the
#'   generated model script, without the `.R` extension. Defaults to
#'   `"biomass_diversity_model_nimble"`.
#'
#' @return Invisibly returns a named list containing:
#' \describe{
#'   \item{model}{Character string containing the generated NIMBLE model code.}
#'   \item{file}{Path to the generated R script.}
#'   \item{response}{Selected response specification.}
#'   \item{biomass_distribution}{Selected biomass observation distribution.}
#'   \item{diversity_distribution}{Selected diversity observation distribution.}
#'   \item{covariates}{Names of the supplied covariates.}
#'   \item{continuous_covariates}{Names of covariates treated as continuous.}
#'   \item{categorical_covariates}{Names of covariates treated as categorical.}
#'   \item{covariate_types}{Named covariate-type specification.}
#'   \item{random_year}{Whether a year random effect was requested.}
#'   \item{random_locality}{Whether a locality random effect was requested.}
#'   \item{random_replicate}{Whether a replicate random effect was requested.}
#'   \item{temporal_effect}{Selected temporal structure.}
#'   \item{temporal_basis}{Original temporal basis matrix.}
#'   \item{temporal_penalty}{Temporal penalty matrix after symmetrization.}
#'   \item{temporal_basis_penalized}{Temporal basis rotated into the retained
#'     penalized eigenbasis, or `NULL` when no cyclic spline is used.}
#'   \item{temporal_penalty_eigenvalues}{Retained positive eigenvalues of the
#'     temporal penalty matrix, or `NULL` when no cyclic spline is used.}
#'   \item{temporal_penalty_eigenvectors}{Eigenvectors associated with the
#'     retained positive eigenvalues, or `NULL` when no cyclic spline is used.}
#'   \item{temporal_penalty_tol}{Tolerance used to retain penalty eigenvalues.}
#'   \item{n_temporal_basis}{Number of retained penalized spline basis
#'     functions.}
#'   \item{temporal_constants}{List of constants required by the penalized
#'     temporal component of the generated model.}
#'   \item{spatial_effect}{Selected spatial structure.}
#'   \item{spatial_shared}{Whether the spatial process is shared.}
#'   \item{spatial_range_prior}{Bounds used for the spatial range prior.}
#'   \item{spatial_sigma_prior}{Scale used for the spatial standard-deviation
#'     prior.}
#' }
#'
#' @details
#' The generated model uses normal observation models for normally distributed
#' biomass and log-normal observation models for positive biomass values.
#' Diversity is modelled using either a Poisson or a negative-binomial
#' likelihood. Continuous covariate effects and intercepts use normal priors.
#' Group-level random effects use zero-centred normal distributions with
#' estimated precisions. Spatial autocorrelation is represented using an
#' exponential distance-based covariance function.
#'
#' When a cyclic cubic spline is requested, the supplied penalty matrix is
#' decomposed as \(S = V D V^\top\). Eigenvectors associated with positive
#' eigenvalues are retained, and the basis is rotated into this eigenbasis.
#' The corresponding spline coefficients receive independent normal priors
#' whose precisions are the product of a response-specific smoothing precision
#' and the retained penalty eigenvalues.
#'
#' This function generates model code but does not compile or fit the model.
#' The caller must supply compatible NIMBLE data, constants, initial values,
#' and MCMC configuration when fitting the model.
#'
#' @examples
#' \dontrun{
#' # Generate a model for biomass and diversity.
#' model_spec <- create_insect_nimble_model(
#'   modelCovs = c("temperature", "habitat"),
#'   response = "both",
#'   biomass_distribution = "lognormal",
#'   diversity_distribution = "negative_binomial",
#'   covariate_types = c(
#'     temperature = "continuous",
#'     habitat = "categorical"
#'   ),
#'   output_dir = tempdir()
#' )
#'
#' # Inspect the generated model and output path.
#' cat(model_spec$model)
#' model_spec$file
#'
#' # Generate a model with a penalized cyclic temporal effect.
#' # temporal_basis and temporal_penalty must be prepared beforehand.
#' cyclic_spec <- create_insect_nimble_model(
#'   modelCovs = "temperature",
#'   response = "diversity",
#'   temporal_effect = "cubic_cyclic",
#'   temporal_basis = temporal_basis,
#'   temporal_penalty = temporal_penalty,
#'   output_dir = tempdir()
#' )
#' }
#'
#' @seealso [nimble::nimbleCode()], [mgcv::cSplineDes()]
#' @export
create_insect_nimble_model <- function(

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

    # --------------------------------------------------------------------------
    # Random effects
    # --------------------------------------------------------------------------

    random_year = FALSE,

    random_locality = FALSE,

    random_replicate = FALSE,

    # --------------------------------------------------------------------------
    # Temporal effect
    # --------------------------------------------------------------------------

    temporal_effect = c(
      "none",
      "linear",
      "cubic_cyclic"
    ),

    temporal_basis = NULL,

    temporal_penalty = NULL,

    temporal_penalty_tol = 1e-8,

    # --------------------------------------------------------------------------
    # Spatial effect
    # --------------------------------------------------------------------------

    spatial_effect = c(
      "none",
      "coordinates",
      "autocorrelation"
    ),

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

    model_name = "biomass_diversity_model_nimble"

) {

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
'
    )


    if (include_biomass) {

      biomass_spatial_terms <- c(
        biomass_spatial_terms,

        "beta_longitude_spatial * longitude[i]",

        "beta_latitude_spatial * latitude[i]"
      )

    }


    if (include_diversity) {

      diversity_spatial_terms <- c(
        diversity_spatial_terms,

        "beta_longitude_spatial * longitude[i]",

        "beta_latitude_spatial * latitude[i]"
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

lambda_spatial_biomass ~ dnorm(
  0,
  0.001
)

lambda_spatial_diversity ~ dnorm(
  0,
  0.001
)
'
      )


      if (include_biomass) {

        biomass_spatial_terms <- c(
          biomass_spatial_terms,

          "lambda_spatial_biomass * spatial_effect[spatial_id[i]]"
        )

      }


      if (include_diversity) {

        diversity_spatial_terms <- c(
          diversity_spatial_terms,

          "lambda_spatial_diversity * spatial_effect[spatial_id[i]]"
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
  # 14. Categorical covariates as random effects
  # ============================================================================

  if (length(categorical_covariates) > 0) {

    for (cov in categorical_covariates) {

      cov_clean <- make.names(cov)


      # ------------------------------------------------------------------------
      # Biomass categorical random effect
      # ------------------------------------------------------------------------

      if (include_biomass) {

        random_effect_code <- c(
          random_effect_code,

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

    if (include_biomass) {

      random_effect_code <- c(
        random_effect_code,

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
      )


      biomass_random_terms <- c(
        biomass_random_terms,

        "year_effect_biomass[year_id[i]]"
      )

    }


    if (include_diversity) {

      random_effect_code <- c(
        random_effect_code,

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

    if (include_biomass) {

      random_effect_code <- c(
        random_effect_code,

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
      )


      biomass_random_terms <- c(
        biomass_random_terms,

        "locality_effect_biomass[locality_id[i]]"
      )

    }


    if (include_diversity) {

      random_effect_code <- c(
        random_effect_code,

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

    if (include_biomass) {

      random_effect_code <- c(
        random_effect_code,

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
      )


      biomass_random_terms <- c(
        biomass_random_terms,

        "replicate_effect_biomass[replicate_id[i]]"
      )

    }


    if (include_diversity) {

      random_effect_code <- c(
        random_effect_code,

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
'
      )


      biomass_temporal_term <- c(
        biomass_temporal_term,

        "beta_time_biomass * temporal[i]"
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
'
      )


      diversity_temporal_term <- c(
        diversity_temporal_term,

        "beta_time_diversity * temporal[i]"
      )

    }

  }


  # ----------------------------------------------------------------------------
  # Penalized cubic cyclic spline
  # ----------------------------------------------------------------------------

  if (temporal_effect == "cubic_cyclic") {

    if (include_biomass) {

      temporal_code <- c(
        temporal_code,

        '
# =============================================================================
# Penalized cubic cyclic spline: biomass
#
# beta_k ~ Normal(
#   0,
#   [tau_temporal * d_k]^-1
# )
#
# where d_k are the positive eigenvalues of the mgcv cyclic-spline
# penalty matrix.
# =============================================================================

for (k in 1:n_temporal_basis) {

  beta_time_biomass[k] ~ dnorm(
    0,
    temporal_tau_biomass *
      temporal_penalty_eigenvalues[k]
  )

}

temporal_tau_biomass ~ dgamma(
  1,
  1
)
'
      )


      biomass_temporal_term <- c(
        biomass_temporal_term,

        paste0(
          "inprod(",
          "beta_time_biomass[1:n_temporal_basis], ",
          "temporal_basis_penalized[",
          "i, 1:n_temporal_basis",
          "])"
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
    temporal_tau_diversity *
      temporal_penalty_eigenvalues[k]
  )

}

temporal_tau_diversity ~ dgamma(
  1,
  1
)
'
      )


      diversity_temporal_term <- c(
        diversity_temporal_term,

        paste0(
          "inprod(",
          "beta_time_diversity[1:n_temporal_basis], ",
          "temporal_basis_penalized[",
          "i, 1:n_temporal_basis",
          "])"
        )
      )

    }

  }


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
  # 23. Fixed-effect priors
  # ============================================================================

  beta_code <- character(0)

  if (n_continuous > 0) {

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

# Temporal effect:
#   ', temporal_effect, '

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
      # Temporal effect
      # ------------------------------------------------------------------------

      temporal_effect = temporal_effect,

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

      spatial_shared = spatial_shared,

      spatial_range_prior = spatial_range_prior,

      spatial_sigma_prior = spatial_sigma_prior

    )
  )

}