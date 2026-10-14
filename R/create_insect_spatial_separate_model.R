
#' Create a NIMBLE Model with Separate Spatial Effects for Insect Biomass and Diversity
#'
#' Constructs a NIMBLE model specification for insect biomass, diversity, or
#' both responses. The model supports continuous and categorical covariates,
#' optional year, locality, and replicate random effects, and optional temporal
#' and spatial effects.
#'
#' Spatial effects can be represented using polynomial coordinate terms or
#' latent spatial autocorrelation processes. For spatial autocorrelation,
#' separate response-specific spatial processes can be fitted, or a shared
#' latent spatial process can be used for both responses. When both responses
#' are included in a shared spatial model, the biomass loading is fixed to
#' one and the diversity loading is estimated to identify the scale of the
#' shared process.
#'
#' @param modelCovs Character vector of covariate names to include in the
#'   model. Covariates are treated as continuous unless `covariate_types`
#'   specifies otherwise.
#'
#' @param response Response to model: `"biomass"`, `"diversity"`, or `"both"`.
#'   Defaults to `"biomass"`.
#'
#' @param biomass_distribution Observation distribution for biomass:
#'   `"normal"` or `"lognormal"`. Defaults to `"normal"`.
#'
#' @param diversity_distribution Observation distribution for diversity:
#'   `"poisson"` or `"negative_binomial"`. Defaults to `"poisson"`.
#'
#' @param covariate_types Optional named character vector specifying whether
#'   each covariate is `"continuous"` or `"categorical"`. Names must correspond
#'   to entries in `modelCovs`. If omitted, all covariates are treated as
#'   continuous.
#'
#' @param random_year Logical; if `TRUE`, includes year-level random effects.
#'   Biomass and diversity have separate year effects and precision parameters.
#'   Defaults to `FALSE`.
#'
#' @param random_locality Logical; if `TRUE`, includes locality-level random
#'   effects. Biomass and diversity have separate locality effects and precision
#'   parameters. Defaults to `FALSE`.
#'
#' @param random_replicate Logical; if `TRUE`, includes replicate-level random
#'   effects. Biomass and diversity have separate replicate effects and
#'   precision parameters. Defaults to `FALSE`.
#'
#' @param temporal_effect Temporal effect to include: `"none"`, `"linear"`,
#'   or `"cubic_cyclic"`. The linear option includes linear and quadratic
#'   temporal terms. The cubic cyclic option uses a supplied spline basis and
#'   penalty matrix. Defaults to `"none"`.
#'
#' @param temporal_basis Numeric matrix of spline basis values for the cubic
#'   cyclic temporal effect. Required when `temporal_effect = "cubic_cyclic"`.
#'   Rows must correspond to observations in the modelling data.
#'
#' @param temporal_penalty Numeric square penalty matrix for the cubic cyclic
#'   spline. Required when `temporal_effect = "cubic_cyclic"` and must have
#'   dimensions equal to the number of columns in `temporal_basis`.
#'
#' @param temporal_penalty_tol Positive, finite numeric scalar used as a
#'   relative threshold for retaining positive eigenvalues of the temporal
#'   penalty matrix. Defaults to `1e-8`.
#'
#' @param temporal_shared Logical; if `TRUE`, uses a shared latent temporal
#'   spline process for both responses. Currently supported only when
#'   `temporal_effect = "cubic_cyclic"`. When both responses are included,
#'   the biomass loading is fixed to one and the diversity loading is
#'   estimated. Defaults to `FALSE`.
#'
#' @param spatial_effect Spatial-effect specification: `"none"`,
#'   `"coordinates"`, or `"autocorrelation"`. The `"coordinates"` option
#'   includes linear and quadratic terms for longitude and latitude.
#'   The `"autocorrelation"` option uses a latent spatial process based on
#'   an exponential distance-decay covariance function. Defaults to
#'   `"autocorrelation"`.
#'
#' @param spatial_shared Logical; if `TRUE` and
#'   `spatial_effect = "autocorrelation"`, uses one latent spatial process
#'   shared by the included responses. If `FALSE`, separate spatial processes
#'   are generated for biomass and diversity. Defaults to `FALSE`.
#'
#' @param spatial_range_prior Numeric vector of length two specifying the
#'   lower and upper bounds for the uniform prior on the spatial range
#'   parameter. Both values must be positive, and the upper bound must
#'   exceed the lower bound. Defaults to `c(1000, 500000)`.
#'
#' @param spatial_sigma_prior Positive, finite numeric scalar specifying the
#'   prior scale for the spatial standard deviation. It determines the
#'   precision used in the half-normal prior for the spatial standard
#'   deviation. Defaults to `1`.
#'
#' @param output_dir Directory in which to write the generated NIMBLE model
#'   file. Defaults to `temp_dir`, which must exist in the calling environment.
#'
#' @param model_name Character string used as the base name of the generated
#'   model file. The output path is constructed as
#'   `file.path(output_dir, paste0(model_name, ".R"))`.
#'   Defaults to `"create_spatial_separate_model"`.
#'
#' @details
#' The function constructs a NIMBLE model as character text, writes the model
#' specification to an R file, and invisibly returns the generated model and
#' associated metadata.
#'
#' **Fixed effects**
#'
#' Continuous covariate effects are estimated separately for biomass and
#' diversity. Each response has its own regression coefficient vector with
#' zero-centred normal priors. Categorical covariates are represented by
#' response-specific random effects, each with its own precision parameter.
#'
#' **Observation models**
#'
#' Biomass can follow a normal or lognormal distribution. Diversity can follow
#' a Poisson or negative-binomial distribution. The model uses separate
#' intercepts and response-specific observation-model parameters.
#'
#' **Random effects**
#'
#' When enabled, year, locality, and replicate effects are generated
#' separately for biomass and diversity. Each response-specific random effect
#' has its own precision parameter.
#'
#' **Spatial coordinate effects**
#'
#' When `spatial_effect = "coordinates"`, the model includes linear and
#' quadratic terms for longitude and latitude. The coordinate coefficients
#' are shared between the included responses in the generated model.
#'
#' **Spatial autocorrelation**
#'
#' When `spatial_effect = "autocorrelation"`, the latent spatial process uses
#' an exponential covariance function:
#'
#' \deqn{
#'   \Sigma_{ij} = \sigma^2
#'   \exp\left(-\frac{d_{ij}}{\rho}\right),
#' }
#'
#' where \eqn{d_{ij}} is the distance between spatial locations,
#' \eqn{\rho} is the spatial range parameter, and \eqn{\sigma} is the
#' spatial standard deviation.
#'
#' The range parameter has a uniform prior between the specified bounds.
#' The spatial standard deviation has a half-normal prior whose precision
#' is determined by `spatial_sigma_prior`.
#'
#' With `spatial_shared = TRUE`, the two responses use the same latent
#' spatial process. If both responses are included, the biomass loading is
#' fixed to one for identifiability and the diversity loading is estimated.
#' If only one response is included, its loading is fixed to one.
#'
#' With `spatial_shared = FALSE`, the model generates separate spatial
#' processes, range parameters, standard deviations, and response-specific
#' loading parameters for the included responses.
#'
#' **Temporal effects**
#'
#' Linear temporal effects include linear and quadratic terms in the variable
#' `temporal`. Cubic cyclic spline effects require a numeric basis matrix and
#' a compatible penalty matrix. The penalty matrix is symmetrised and
#' decomposed into eigenvalues and eigenvectors. Eigenvalues above the
#' relative threshold defined by `temporal_penalty_tol` are retained, and
#' the spline basis is rotated into the retained penalized eigenbasis.
#'
#' When `temporal_shared = TRUE`, the cubic cyclic spline coefficients and
#' smoothing precision are shared. For a two-response model, the biomass
#' loading is fixed to one and the diversity loading is estimated.
#' When `temporal_shared = FALSE`, response-specific spline coefficients
#' and smoothing precision parameters are generated.
#'
#' The generated model expects all required NIMBLE data, dimensions, indices,
#' and constants to be supplied when the model is instantiated. These include
#' `X` for continuous covariates; `longitude` and `latitude` for coordinate
#' effects; `distance_matrix`, `spatial_zero`, and spatial indices for spatial
#' autocorrelation; and `year_id`, `locality_id`, and `replicate_id` for
#' the corresponding random effects when enabled.
#'
#' For cubic cyclic temporal effects, the returned `temporal_constants` list
#' contains the penalized basis matrix, retained penalty eigenvalues, and
#' number of retained basis functions. These should be supplied as NIMBLE
#' constants.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{model}{Generated NIMBLE model specification as character text.}
#'   \item{file}{Path to the generated model file.}
#'   \item{response}{Selected response or responses.}
#'   \item{biomass_distribution, diversity_distribution}{Selected observation
#'     distributions.}
#'   \item{covariates}{All specified covariates.}
#'   \item{continuous_covariates, categorical_covariates}{Covariates grouped
#'     by type.}
#'   \item{covariate_types}{Named covariate-type specification.}
#'   \item{random_year, random_locality, random_replicate}{Random-effect
#'     inclusion flags.}
#'   \item{effect_structure}{Selected model effect structure.}
#'   \item{temporal_effect, temporal_shared}{Temporal model configuration.}
#'   \item{temporal_basis, temporal_penalty}{Input spline basis and penalty.}
#'   \item{temporal_basis_penalized}{Basis rotated into the penalized
#'     eigenbasis.}
#'   \item{temporal_penalty_eigenvalues, temporal_penalty_eigenvectors}{Retained
#'     penalty eigenvalues and eigenvectors.}
#'   \item{temporal_penalty_tol, n_temporal_basis}{Penalty threshold and number
#'     of retained spline basis functions.}
#'   \item{temporal_constants}{Constants required for the cyclic spline.}
#'   \item{spatial_effect, spatial_shared}{Spatial model configuration.}
#'   \item{spatial_range_prior, spatial_sigma_prior}{Spatial prior settings.}
#' }
#'
#' @examples
#' \dontrun{
#' # Biomass model with separate spatial autocorrelation
#' biomass_model <- create_insect_spatial_separate_model(
#'   modelCovs = c("temperature", "vegetation"),
#'   response = "biomass",
#'   covariate_types = c(
#'     temperature = "continuous",
#'     vegetation = "categorical"
#'   ),
#'   spatial_effect = "autocorrelation",
#'   spatial_shared = FALSE,
#'   output_dir = tempdir(),
#'   model_name = "insect_biomass_spatial"
#' )
#'
#' # Joint model with a shared latent spatial process
#' shared_model <- create_insect_spatial_separate_model(
#'   modelCovs = c("temperature"),
#'   response = "both",
#'   spatial_effect = "autocorrelation",
#'   spatial_shared = TRUE,
#'   output_dir = tempdir(),
#'   model_name = "insect_shared_spatial"
#' )
#'
#' # Joint model with separate spatial processes and cyclic temporal effects
#' spatiotemporal_model <- create_insect_spatial_separate_model(
#'   modelCovs = c("temperature"),
#'   response = "both",
#'   spatial_effect = "autocorrelation",
#'   spatial_shared = FALSE,
#'   temporal_effect = "cubic_cyclic",
#'   temporal_basis = spline_basis,
#'   temporal_penalty = spline_penalty,
#'   temporal_shared = TRUE,
#'   output_dir = tempdir(),
#'   model_name = "insect_spatiotemporal"
#' )
#' }
#'
#' @seealso [nimble::nimbleCode()]
#'
#' @export


create_insect_spatial_separate_model <- function(

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

  spatial_shared = FALSE,

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

  model_name = "create_spatial_separate_model"

) {

  # Fixed structure for this standalone model builder: Spatial-only model with separate response-specific spatial processes.
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

# Loadings are added conditionally below so one reference loading can be fixed.
'
      )


if (include_biomass) {

  # Reference loading fixed to 1 for identifiability.
  biomass_spatial_terms <- c(
    biomass_spatial_terms,
    "spatial_effect[spatial_id[i]]"
  )

}

if (include_diversity) {

  if (include_biomass) {
    spatial_code <- c(
      spatial_code,
      '
# Estimate only the diversity loading; biomass loading is fixed to 1.
lambda_spatial_diversity ~ dnorm(0, 1)
'
    )
    diversity_spatial_terms <- c(
      diversity_spatial_terms,
      "lambda_spatial_diversity * spatial_effect[spatial_id[i]]"
    )
  } else {
    # Single-response model: loading fixed to 1.
    diversity_spatial_terms <- c(
      diversity_spatial_terms,
      "spatial_effect[spatial_id[i]]"
    )
  }

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
