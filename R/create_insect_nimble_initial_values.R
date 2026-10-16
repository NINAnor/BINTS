#' Create Initial Values for the NIMBLE Insect Model
#'
#' Generate conservative initial values for parameters in a NIMBLE model of
#' insect biomass, diversity, or both responses. Initial values are constructed
#' from the observed response data and the model configuration so that parameter
#' names and vector dimensions correspond to the selected model components.
#'
#' The function supports normal or lognormal biomass models, Poisson or negative
#' binomial diversity models, continuous and categorical covariates, optional
#' year, locality, and replicate random effects, linear or cubic cyclic temporal
#' effects, and coordinate-based or autocorrelated spatial effects.
#'
#' @param model_structure Optional explicit model structure: `baseline`,
#'   `temporal_shared`, `temporal_separate`, `spatial_shared`,
#'   `spatial_separate`, `spatiotemporal_shared`, or `spatiotemporal_separate`.
#'   If `NULL` (default), the structure is inferred from the temporal and spatial
#'   arguments, like the model-building dispatcher.
#' @param modelling_data A data frame containing the observations, response
#'   variables, covariates, and any grouping or spatial coordinate columns
#'   required by the selected model configuration.
#' @param modelCovs Character vector of covariate column names included in the
#'   model. If `NULL` or omitted, no covariates are included.
#' @param covariate_types Optional named character vector specifying the type of
#'   each covariate in `modelCovs`. Supported values are `"continuous"` and
#'   `"categorical"`. Names must cover all specified covariates. If omitted,
#'   all covariates are treated as continuous.
#' @param response Character string specifying the response or responses for
#'   which initial values are required: `"biomass"`, `"diversity"`, or `"both"`.
#' @param biomass_col Name of the column containing the biomass response.
#'   Defaults to `"sum_wet_weight"`. Required when `response` includes biomass.
#' @param diversity_col Name of the column containing the diversity response.
#'   Required when `response` includes diversity.
#' @param biomass_distribution Biomass response distribution: `"normal"` or
#'   `"lognormal"`. Defaults to `"normal"`.
#' @param diversity_distribution Diversity response distribution: `"poisson"`
#'   or `"negative_binomial"`. Defaults to `"poisson"`.
#' @param random_year Logical; if `TRUE`, initialise one shared year-level
#'   random effect used by each included response.
#' @param random_locality Logical; if `TRUE`, initialise one shared locality-
#'   level random effect used by each included response.
#' @param random_replicate Logical; if `TRUE`, initialise one shared replicate-
#'   level random effect used by each included response. Replicate numbers are
#'   extracted from the suffix following the last hyphen in the replicate
#'   identifier and indexed within locality-year groups.
#' @param year_col Name of the column identifying years when `random_year =
#'   TRUE`. Defaults to `"year"`.
#' @param locality_col Name of the column identifying localities when
#'   `random_locality = TRUE`. Defaults to `"locality"`.
#' @param replicate_col Name of the column containing replicate identifiers
#'   when `random_replicate = TRUE`. Defaults to `"sampling_name"`.
#' @param temporal_data Numeric vector containing the temporal covariate, with
#'   one value per row of `modelling_data`. Required when `temporal_effect` is
#'   `"linear"` or `"cubic_cyclic"`.
#' @param temporal_effect Temporal specification: `"none"`, `"linear"`, or
#'   `"cubic_cyclic"`. When `model_structure = NULL`, this determines whether
#'   the temporal component is included. Defaults to `"none"`. If omitted, the selected `model_structure` determines
#'   whether a temporal effect is included and defaults to `"cubic_cyclic"`.
#' @param temporal_shared Logical or `NULL`. Normally inferred from
#'   `model_structure`; shared temporal coefficients are supported only for
#'   cubic cyclic effects. If supplied, the value must agree with the selected
#'   model structure.
#' @param temporal_basis_penalized Numeric matrix containing the cyclic spline
#'   basis projected onto the positive-eigenvalue subspace of the spline penalty
#'   matrix. Required when `temporal_effect = "cubic_cyclic"`. It must have one
#'   row per observation and at least one column.
#' @param spatial_effect Optional spatial specification: `"coordinates"` or
#'   `"autocorrelation"`. If omitted, spatial model structures default to
#'   `"autocorrelation"` and other structures use no spatial effect.
#' @param spatial_shared Optional logical override for the spatial process. The
#'   selected structure sets this automatically; it is exposed for compatibility
#'   and must agree with `model_structure`.
#' @param longitude_col Name of the numeric longitude column. Defaults to
#'   `"longitude"`.
#' @param latitude_col Name of the numeric latitude column. Defaults to
#'   `"latitude"`.
#' @param spatial_range_prior Numeric vector of length two giving the lower and
#'   upper bounds of the spatial range prior used by the model. Both values must
#'   be positive and finite, and the upper bound must exceed the lower bound.
#'   Defaults to `c(1000, 500000)`.
#' @param spatial_sigma_prior Positive finite numeric value specifying the
#'   spatial standard-deviation prior parameter used by the model. Defaults to
#'   `1`.
#' @param seed Optional integer or numeric value passed to [set.seed()] for
#'   reproducible initial-value generation. If `NULL`, the random-number
#'   generator is not explicitly seeded.
#'
#' @details
#' Initial values are returned as a named list whose elements correspond to
#' parameters in the generated NIMBLE model.
#'
#' For biomass, the intercept is initialised to the observed mean for a normal
#' model and to the mean log biomass for a lognormal model. The biomass
#' precision is initialised using the inverse squared sample standard deviation
#' on the relevant response scale, with a fallback value when the standard
#' deviation is non-finite or non-positive.
#'
#' For Poisson and negative binomial diversity models, the intercept is
#' initialised on the log scale using the observed mean. A small positive
#' fallback is used when the mean is zero or negative. For negative binomial
#' models, the dispersion parameter is initialised using a method-of-moments
#' estimate when the observed variance exceeds the mean; otherwise, a fallback
#' value is used. The resulting estimate is constrained to the interval
#' [0.1, 1000].
#'
#' The selected model structure determines whether fixed and ordinary random
#' effects are shared or response-specific. Initial values are zero for regression
#' coefficients and latent random effects, with precision parameters set to one.
#'
#' For a linear temporal effect, each included response receives a temporal
#' regression coefficient initialised to zero. For a cubic cyclic temporal
#' effect, each response receives its own zero-initialised coefficient vector
#' and precision by default. When `temporal_shared = TRUE`, a single shared
#' coefficient vector (`beta_time_shared`) and precision
#' (`temporal_tau_shared`) are initialised instead. The supplied basis must be
#' consistent with the penalised basis used by the model-construction function.
#'
#' For a coordinate-based spatial effect, longitude and latitude coefficients
#' are initialised to zero. For spatial autocorrelation, latent spatial effects
#' are initialised to zero, range parameters to the midpoint of
#' `spatial_range_prior`, and spatial standard deviations to the larger of 0.5
#' and `spatial_sigma_prior`. Separate spatial structures receive response-
#' specific fields. The spatiotemporal shared builder also initialises the
#' diversity loading for its shared spatial field.
#'
#' The response columns are checked for numeric, finite values. Lognormal
#' biomass requires strictly positive observations. Poisson and negative
#' binomial diversity require non-negative integer observations. Spatial
#' coordinates must be numeric, finite, non-missing, and within valid
#' longitude and latitude ranges.
#'
#' The selected configuration must agree with the model code and constants
#' generated for fitting. In particular, temporal spline dimensions and spatial
#' process dimensions must match those used by the other model-construction
#' helpers.
#'
#' @return A named list of initial parameter values suitable for use as the
#'   `inits` argument when configuring a NIMBLE model or MCMC sampler. The
#'   returned elements depend on the selected responses and model components.
#'
#' @seealso
#' [create_insect_spatial_shared_model()], [create_insect_nimble_data()],
#' [create_insect_nimble_constants()]
#'
#' @export


create_insect_nimble_initial_values <- function(

  modelling_data,
  modelCovs,

  model_structure = NULL,
  covariate_types = NULL,
  share_fixed_effects = FALSE,
    share_hyperparameters = FALSE,

  # --------------------------------------------------------------------------
  # Response specification
  # --------------------------------------------------------------------------

  response = "biomass",

  biomass_col = "sum_wet_weight",

  diversity_col = NULL,

  # --------------------------------------------------------------------------
  # Response distributions
  # --------------------------------------------------------------------------

  biomass_distribution = "normal",

  diversity_distribution = "poisson",

  # --------------------------------------------------------------------------
  # Random effects
  # --------------------------------------------------------------------------

  random_year = FALSE,

  random_locality = FALSE,

  random_replicate = FALSE,

  # --------------------------------------------------------------------------
  # Column names
  # --------------------------------------------------------------------------

  year_col = "year",

  locality_col = "locality",

  replicate_col = "sampling_name",

  # --------------------------------------------------------------------------
  # Temporal effect
  # --------------------------------------------------------------------------

  temporal_data = NULL,

  temporal_effect = NULL,

  temporal_shared = NULL,

  # --------------------------------------------------------------------------
  # Temporal spline basis
  #
  # For cubic_cyclic this should be the already transformed basis:
  #
  #     temporal_basis_penalized
  #
  # generated by create_insect_nimble_model().
  #
  # --------------------------------------------------------------------------

  temporal_basis_penalized = NULL,

  # --------------------------------------------------------------------------
  # Spatial effect
  # --------------------------------------------------------------------------

  spatial_effect = NULL,

  spatial_shared = NULL,

  longitude_col = "longitude",

  latitude_col = "latitude",

  # --------------------------------------------------------------------------
  # Spatial parameter priors
  #
  # These should match the priors used in
  # create_insect_nimble_model().
  # --------------------------------------------------------------------------

  spatial_range_prior = c(
    1000,
    500000
  ),

  spatial_sigma_prior = 1,

  # --------------------------------------------------------------------------
  # Reproducibility
  # --------------------------------------------------------------------------

  seed = NULL

) {


  # ==========================================================================
  # 1. Set seed
  # ==========================================================================

  if (!is.null(seed)) {

    set.seed(
      seed
    )

  }


  # ==========================================================================
  # ==========================================================================
  # 2. Resolve the seven-model structure
  # ==========================================================================
  # Infer the structure from the same arguments used by the model dispatcher
  # unless the caller explicitly supplies model_structure.

  available_structures <- c(
    "baseline", "temporal_shared", "temporal_separate",
    "spatial_shared", "spatial_separate",
    "spatiotemporal_shared", "spatiotemporal_separate"
  )
  structure_was_explicit <- !is.null(model_structure)

  if (structure_was_explicit) {
    model_structure <- match.arg(model_structure, choices = available_structures)
  } else {
    if (is.null(temporal_effect)) temporal_effect <- "none"
    if (is.null(spatial_effect)) spatial_effect <- "none"
    if (is.null(temporal_shared)) temporal_shared <- TRUE
    if (is.null(spatial_shared)) spatial_shared <- TRUE

    if (!temporal_effect %in% c("none", "linear", "cubic_cyclic")) {
      stop("'temporal_effect' must be 'none', 'linear', or 'cubic_cyclic'.",
           call. = FALSE)
    }
    if (!spatial_effect %in% c("none", "coordinates", "autocorrelation")) {
      stop("'spatial_effect' must be 'none', 'coordinates', or 'autocorrelation'.",
           call. = FALSE)
    }
    if (!is.logical(temporal_shared) || length(temporal_shared) != 1L ||
        is.na(temporal_shared)) {
      stop("'temporal_shared' must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(spatial_shared) || length(spatial_shared) != 1L ||
        is.na(spatial_shared)) {
      stop("'spatial_shared' must be TRUE or FALSE.", call. = FALSE)
    }

    has_temporal <- temporal_effect != "none"
    has_spatial <- spatial_effect != "none"

    if (!has_temporal && !has_spatial) {
      model_structure <- "baseline"
    } else if (has_temporal && !has_spatial) {
      model_structure <- if (temporal_shared) "temporal_shared" else "temporal_separate"
    } else if (!has_temporal && has_spatial) {
      model_structure <- if (spatial_shared) "spatial_shared" else "spatial_separate"
    } else {
      if (temporal_shared != spatial_shared) {
        stop(
          "Spatiotemporal models require temporal_shared and ",
          "spatial_shared to have the same value.",
          call. = FALSE
        )
      }
      model_structure <- if (temporal_shared) {
        "spatiotemporal_shared"
      } else {
        "spatiotemporal_separate"
      }
    }
  }

  structure_defaults <- switch(
    model_structure,
    baseline = list(temporal = "none", spatial = "none",
                    temporal_shared = FALSE, spatial_shared = TRUE),
    temporal_shared = list(temporal = "cubic_cyclic", spatial = "none",
                           temporal_shared = TRUE, spatial_shared = TRUE),
    temporal_separate = list(temporal = "cubic_cyclic", spatial = "none",
                             temporal_shared = FALSE, spatial_shared = TRUE),
    spatial_shared = list(temporal = "none", spatial = "autocorrelation",
                          temporal_shared = FALSE, spatial_shared = TRUE),
    spatial_separate = list(temporal = "none", spatial = "autocorrelation",
                            temporal_shared = FALSE, spatial_shared = FALSE),
    spatiotemporal_shared = list(temporal = "cubic_cyclic",
                                 spatial = "autocorrelation",
                                 temporal_shared = TRUE, spatial_shared = TRUE),
    spatiotemporal_separate = list(temporal = "cubic_cyclic",
                                   spatial = "autocorrelation",
                                   temporal_shared = FALSE, spatial_shared = FALSE)
  )

  if (is.null(temporal_effect)) temporal_effect <- structure_defaults$temporal
  if (is.null(spatial_effect)) spatial_effect <- structure_defaults$spatial
  if (is.null(temporal_shared)) temporal_shared <- structure_defaults$temporal_shared
  if (is.null(spatial_shared)) spatial_shared <- structure_defaults$spatial_shared

  # Ignore sharing flags for effects that are absent, as the model dispatcher
  # does. For present effects, the flags must match the selected structure.
  if (temporal_effect == "none") temporal_shared <- structure_defaults$temporal_shared
  if (spatial_effect == "none") spatial_shared <- structure_defaults$spatial_shared

  if (!identical(isTRUE(temporal_shared), structure_defaults$temporal_shared)) {
    stop("'temporal_shared' must agree with model_structure = '",
         model_structure, "'.", call. = FALSE)
  }
  if (!identical(isTRUE(spatial_shared), structure_defaults$spatial_shared)) {
    stop("'spatial_shared' must agree with model_structure = '",
         model_structure, "'.", call. = FALSE)
  }
  if (model_structure %in% c("baseline", "spatial_shared", "spatial_separate") &&
      temporal_effect != "none") {
    stop("model_structure = '", model_structure,
         "' does not include a temporal effect.", call. = FALSE)
  }
  if (model_structure %in% c("baseline", "temporal_shared", "temporal_separate") &&
      spatial_effect != "none") {
    stop("model_structure = '", model_structure,
         "' does not include a spatial effect.", call. = FALSE)
  }
  if (model_structure %in% c("temporal_shared", "temporal_separate",
                             "spatiotemporal_shared", "spatiotemporal_separate") &&
      temporal_effect == "none") {
    stop("model_structure = '", model_structure,
         "' requires a temporal effect.", call. = FALSE)
  }
  if (model_structure %in% c("spatial_shared", "spatial_separate",
                             "spatiotemporal_shared", "spatiotemporal_separate") &&
      spatial_effect == "none") {
    stop("model_structure = '", model_structure,
         "' requires a spatial effect.", call. = FALSE)
  }
  if (model_structure %in% c("spatial_separate", "spatiotemporal_separate") &&
      spatial_effect == "coordinates") {
    stop("Separate spatial builders support spatial autocorrelation, not the shared coordinate trend.",
         call. = FALSE)
  }

  # 3. Validate response
  # ==========================================================================

  response <- match.arg(
    response,
    choices = c(
      "biomass",
      "diversity",
      "both"
    )
  )


  # ==========================================================================
  # 3. Validate response distributions
  # ==========================================================================

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


  # ==========================================================================
  # 4. Validate temporal effect
  # ==========================================================================

  temporal_effect <- match.arg(
    temporal_effect,
    choices = c(
      "none",
      "linear",
      "cubic_cyclic"
    )
  )

  if (
    !is.logical(temporal_shared) ||
    length(temporal_shared) != 1 ||
    is.na(temporal_shared)
  ) {
    stop("'temporal_shared' must be TRUE or FALSE.", call. = FALSE)
  }

  # if (isTRUE(temporal_shared) && temporal_effect != "cubic_cyclic") {
  #   stop(
  #     "'temporal_shared = TRUE' is supported only when ",
  #     "temporal_effect = 'cubic_cyclic'.",
  #     call. = FALSE
  #   )
  # }


  # ==========================================================================
  # 5. Validate spatial effect
  # ==========================================================================

  spatial_effect <- match.arg(
    spatial_effect,
    choices = c(
      "none",
      "coordinates",
      "autocorrelation"
    )
  )


  # ==========================================================================
  # 6. Validate spatial_shared
  # ==========================================================================

  if (
    !is.logical(spatial_shared) ||
    length(spatial_shared) != 1 ||
    is.na(spatial_shared)
  ) {
    stop("'spatial_shared' must be TRUE or FALSE.", call. = FALSE)
  }



  # ==========================================================================
  # 7. Validate spatial parameter priors
  # ==========================================================================

  if (
    !is.numeric(spatial_range_prior) ||
    length(spatial_range_prior) != 2 ||
    any(!is.finite(spatial_range_prior)) ||
    spatial_range_prior[1] <= 0 ||
    spatial_range_prior[2] <= spatial_range_prior[1]
  ) {

    stop(
      "'spatial_range_prior' must contain two positive finite values ",
      "with the lower bound smaller than the upper bound.",
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
      "'spatial_sigma_prior' must be a single positive finite value.",
      call. = FALSE
    )

  }


  # ==========================================================================
  # 8. Determine included responses
  # ==========================================================================

  include_biomass <- response %in% c(
    "biomass",
    "both"
  )

  include_diversity <- response %in% c(
    "diversity",
    "both"
  )


  # ==========================================================================
  # 9. Validate modelling data
  # ==========================================================================

  if (!is.data.frame(modelling_data)) {

    stop(
      "modelling_data must be a data frame.",
      call. = FALSE
    )

  }

  N <- nrow(
    modelling_data
  )

  if (N == 0) {

    stop(
      "modelling_data contains no observations.",
      call. = FALSE
    )

  }


  # ==========================================================================
  # 10. Validate response columns
  # ==========================================================================

  if (include_biomass) {

    if (is.null(biomass_col)) {

      stop(
        "biomass_col must be supplied when ",
        "response includes biomass.",
        call. = FALSE
      )

    }

    if (!biomass_col %in% names(modelling_data)) {

      stop(
        "Biomass column not found in modelling_data: ",
        biomass_col,
        call. = FALSE
      )

    }

  }


  if (include_diversity) {

    if (is.null(diversity_col)) {

      stop(
        "diversity_col must be supplied when ",
        "response includes diversity.",
        call. = FALSE
      )

    }

    if (!diversity_col %in% names(modelling_data)) {

      stop(
        "Diversity column not found in modelling_data: ",
        diversity_col,
        call. = FALSE
      )

    }

  }


  # ==========================================================================
  # 11. Validate model covariates
  # ==========================================================================

  if (missing(modelCovs) || is.null(modelCovs)) {

    modelCovs <- character(0)

  }

  if (!is.character(modelCovs)) {

    stop(
      "modelCovs must be a character vector.",
      call. = FALSE
    )

  }

  if (anyDuplicated(modelCovs)) {

    stop(
      "modelCovs contains duplicated covariate names.",
      call. = FALSE
    )

  }


  # ==========================================================================
  # 12. Validate covariate types
  # ==========================================================================

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

      if (!is.character(covariate_types)) {

        stop(
          "covariate_types must be a character vector.",
          call. = FALSE
        )

      }

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
          "'continuous' or 'categorical'. ",
          "Invalid type(s): ",
          paste(
            invalid_types,
            collapse = ", "
          ),
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


  # ==========================================================================
  # 13. Separate covariates
  # ==========================================================================

  continuous_covariates <- modelCovs[
    covariate_types[modelCovs] == "continuous"
  ]

  categorical_covariates <- modelCovs[
    covariate_types[modelCovs] == "categorical"
  ]


  # ==========================================================================
  # 14. Validate spatial coordinates
  # ==========================================================================

  if (
    spatial_effect %in%
    c(
      "coordinates",
      "autocorrelation"
    )
  ) {

    required_spatial_cols <- c(
      longitude_col,
      latitude_col
    )

    missing_spatial_cols <- setdiff(
      required_spatial_cols,
      names(modelling_data)
    )

    if (length(missing_spatial_cols) > 0) {

      stop(
        "The following spatial coordinate columns are missing from ",
        "modelling_data: ",
        paste(
          missing_spatial_cols,
          collapse = ", "
        ),
        call. = FALSE
      )

    }

    longitude <- modelling_data[[longitude_col]]
    latitude <- modelling_data[[latitude_col]]

    if (!is.numeric(longitude)) {

      stop(
        "'",
        longitude_col,
        "' must be numeric.",
        call. = FALSE
      )

    }

    if (!is.numeric(latitude)) {

      stop(
        "'",
        latitude_col,
        "' must be numeric.",
        call. = FALSE
      )

    }

    if (
      anyNA(longitude) ||
      anyNA(latitude)
    ) {

      stop(
        "Missing longitude/latitude values are not allowed when ",
        "a spatial effect is used.",
        call. = FALSE
      )

    }

    if (
      any(!is.finite(longitude)) ||
      any(!is.finite(latitude))
    ) {

      stop(
        "Longitude and latitude must contain only finite values.",
        call. = FALSE
      )

    }

    if (
      any(
        longitude < -180 |
        longitude > 180
      )
    ) {

      stop(
        "Longitude values must be between -180 and 180 degrees.",
        call. = FALSE
      )

    }

    if (
      any(
        latitude < -90 |
        latitude > 90
      )
    ) {

      stop(
        "Latitude values must be between -90 and 90 degrees.",
        call. = FALSE
      )

    }

  }


  # ==========================================================================
  # 15. Initialise output
  # ==========================================================================

  initial_values <- list()


  # ==========================================================================
  # 16. Biomass initial values
  # ==========================================================================

  if (include_biomass) {

    biomass <- modelling_data[[biomass_col]]

    # ----------------------------------------------------------------------
    # Validate numeric response
    # ----------------------------------------------------------------------

    if (!is.numeric(biomass)) {

      stop(
        "Biomass response must be numeric.",
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Validate finite values
    # ----------------------------------------------------------------------

    if (anyNA(biomass)) {

      stop(
        "Missing values found in biomass response.",
        call. = FALSE
      )

    }

    if (any(!is.finite(biomass))) {

      stop(
        "Non-finite values found in biomass response.",
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Log-normal biomass requires positive observations
    # ----------------------------------------------------------------------

    if (biomass_distribution == "lognormal") {

      if (any(biomass <= 0)) {

        stop(
          "Log-normal biomass requires strictly positive observations. ",
          "Found ",
          sum(biomass <= 0),
          " observation(s) <= 0.",
          call. = FALSE
        )

      }

    }

    # ----------------------------------------------------------------------
    # Biomass intercept
    # ----------------------------------------------------------------------

    if (biomass_distribution == "normal") {

      initial_values[[
        "alpha_biomass"
      ]] <- mean(
        biomass
      )

    } else {

      # For a log-normal model, alpha is on the log scale.

      initial_values[[
        "alpha_biomass"
      ]] <- mean(
        log(biomass)
      )

    }

    # ----------------------------------------------------------------------
    # Biomass residual precision
    # ----------------------------------------------------------------------

    if (biomass_distribution == "normal") {

      biomass_sd <- stats::sd(
        biomass
      )

    } else {

      biomass_sd <- stats::sd(
        log(biomass)
      )

    }

    if (
      !is.finite(biomass_sd) ||
      biomass_sd <= 0
    ) {

      biomass_sd <- 1

    }

    initial_values[[
      "tau_biomass"
    ]] <- 1 / biomass_sd^2

  }


  # ==========================================================================
  # 17. Diversity initial values
  # ==========================================================================

  if (include_diversity) {

    diversity <- modelling_data[[diversity_col]]

    # ----------------------------------------------------------------------
    # Validate numeric response
    # ----------------------------------------------------------------------

    if (!is.numeric(diversity)) {

      stop(
        "Diversity response must be numeric.",
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Validate missing and non-finite values
    # ----------------------------------------------------------------------

    if (anyNA(diversity)) {

      stop(
        "Missing values found in diversity response.",
        call. = FALSE
      )

    }

    if (any(!is.finite(diversity))) {

      stop(
        "Non-finite values found in diversity response.",
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Validate count response
    # ----------------------------------------------------------------------

    if (
      diversity_distribution %in%
      c(
        "poisson",
        "negative_binomial"
      )
    ) {

      if (any(diversity < 0)) {

        stop(
          diversity_distribution,
          " diversity requires non-negative observations.",
          call. = FALSE
        )

      }

      if (any(diversity != floor(diversity))) {

        stop(
          diversity_distribution,
          " diversity requires integer-valued observations.",
          call. = FALSE
        )

      }

    }

    # ----------------------------------------------------------------------
    # Diversity intercept
    # ----------------------------------------------------------------------

    if (
      diversity_distribution %in%
      c(
        "poisson",
        "negative_binomial"
      )
    ) {

      mean_diversity <- mean(
        diversity
      )

      # Avoid log(0)

      if (mean_diversity <= 0) {

        mean_diversity <- 0.1

      }

      initial_values[[
        "alpha_diversity"
      ]] <- log(
        mean_diversity
      )

    }

    # Continuous-covariate coefficients are shared across responses and
    # are initialised once as `beta` in the biomass section above.

    # ----------------------------------------------------------------------
    # Negative-binomial dispersion
    # ----------------------------------------------------------------------

    if (
      diversity_distribution ==
      "negative_binomial"
    ) {

      # ------------------------------------------------------------------
      # Method-of-moments initial value:
      #
      #   Var(Y) = mu + mu^2 / size
      #
      # therefore:
      #
      #   size = mu^2 / (Var(Y) - mu)
      # ------------------------------------------------------------------

      mean_diversity <- mean(
        diversity
      )

      variance_diversity <- stats::var(
        diversity
      )

      if (
        is.finite(mean_diversity) &&
        is.finite(variance_diversity) &&
        variance_diversity > mean_diversity &&
        mean_diversity > 0
      ) {

        size_diversity <- (
          mean_diversity^2 /
            (
              variance_diversity -
                mean_diversity
            )
        )

      } else {

        size_diversity <- 10

      }

      # ------------------------------------------------------------------
      # Keep initial value numerically reasonable
      # ------------------------------------------------------------------

      size_diversity <- max(
        0.1,
        min(
          size_diversity,
          1000
        )
      )

      initial_values[[
        "size_diversity"
      ]] <- size_diversity

    }

  }


  # ==========================================================================
  # Shared continuous-covariate coefficients
  # ==========================================================================

  if (length(continuous_covariates) > 0) {
    initial_values[["beta"]] <- rep(
      0,
      length(continuous_covariates)
    )
  }


  # ==========================================================================
  # 18. Shared categorical covariate effects
  # ==========================================================================

  if (length(categorical_covariates) > 0) {
    for (cov in categorical_covariates) {
      if (!cov %in% names(modelling_data)) {
        stop(
          "Categorical covariate not found in modelling_data: ",
          cov,
          call. = FALSE
        )
      }

      cov_values <- modelling_data[[cov]]
      if (anyNA(cov_values)) {
        stop(
          "Missing values found in categorical covariate: ",
          cov,
          ".",
          call. = FALSE
        )
      }

      n_levels <- nlevels(factor(cov_values))
      cov_clean <- make.names(cov)
      initial_values[[paste0(cov_clean, "_effect")]] <- rep(0, n_levels)
      initial_values[[paste0("tau_", cov_clean)]] <- 1
    }
  }


  # ==========================================================================
  # 19. Shared year random effect
  # ==========================================================================

  if (random_year) {
    if (!year_col %in% names(modelling_data)) {
      stop(
        "Year column not found in modelling_data: ",
        year_col,
        call. = FALSE
      )
    }
    year_values <- modelling_data[[year_col]]
    if (anyNA(year_values)) {
      stop("Missing values found in year column.", call. = FALSE)
    }

    initial_values[["year_effect"]] <- rep(
      0,
      dplyr::n_distinct(year_values)
    )
    initial_values[["tau_year"]] <- 1
  }


  # ==========================================================================
  # 20. Shared locality random effect
  # ==========================================================================

  if (random_locality) {
    if (!locality_col %in% names(modelling_data)) {
      stop(
        "Locality column not found in modelling_data: ",
        locality_col,
        call. = FALSE
      )
    }
    locality_values <- modelling_data[[locality_col]]
    if (anyNA(locality_values)) {
      stop("Missing values found in locality column.", call. = FALSE)
    }

    initial_values[["locality_effect"]] <- rep(
      0,
      dplyr::n_distinct(locality_values)
    )
    initial_values[["tau_locality"]] <- 1
  }


  # ==========================================================================
  # 21. Replicate random effect
  #
  # Replicate numbers are extracted from replicate_col and then indexed within
  # locality × year, matching create_insect_nimble_constants().
  #
  # ==========================================================================

  if (random_replicate) {

    # ----------------------------------------------------------------------
    # Check replicate column
    # ----------------------------------------------------------------------

    if (is.null(replicate_col)) {

      stop(
        "replicate_col must be supplied when ",
        "random_replicate = TRUE.",
        call. = FALSE
      )

    }

    if (!replicate_col %in% names(modelling_data)) {

      stop(
        "Replicate column not found in modelling_data: ",
        replicate_col,
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Check grouping variables
    # ----------------------------------------------------------------------

    required_replicate_vars <- c(
      "locality",
      "year"
    )

    missing_replicate_vars <- setdiff(
      required_replicate_vars,
      names(modelling_data)
    )

    if (length(missing_replicate_vars) > 0) {

      stop(
        "The following columns are required to construct ",
        "replicate IDs: ",
        paste(
          missing_replicate_vars,
          collapse = ", "
        ),
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Extract replicate number
    # ----------------------------------------------------------------------

    replicate_number <- suppressWarnings(
      as.numeric(
        sub(
          ".*-",
          "",
          as.character(
            modelling_data[[replicate_col]]
          )
        )
      )
    )

    # ----------------------------------------------------------------------
    # Validate extracted replicate numbers
    # ----------------------------------------------------------------------

    if (anyNA(replicate_number)) {

      bad_values <- unique(
        modelling_data[[replicate_col]][
          is.na(replicate_number)
        ]
      )

      stop(
        "Could not extract numeric replicate numbers from '",
        replicate_col,
        "'. Problematic value(s): ",
        paste(
          bad_values,
          collapse = ", "
        ),
        call. = FALSE
      )

    }

    # ----------------------------------------------------------------------
    # Construct replicate IDs within locality × year
    # ----------------------------------------------------------------------

    replicates <- modelling_data |>
      dplyr::transmute(
        locality = .data$locality,
        year = .data$year,
        replicate_number = replicate_number
      ) |>
      dplyr::group_by(
        locality,
        year
      ) |>
      dplyr::mutate(
        sampling_id = as.integer(
          factor(
            replicate_number
          )
        )
      ) |>
      dplyr::ungroup()

    # ----------------------------------------------------------------------
    # Validate replicate IDs
    # ----------------------------------------------------------------------

    if (anyNA(replicates$sampling_id)) {

      stop(
        "Missing values found in replicate IDs.",
        call. = FALSE
      )

    }

    n_replicate <- max(
      replicates$sampling_id
    )

    # ----------------------------------------------------------------------
    # Shared replicate effect
    # ----------------------------------------------------------------------

    initial_values[["replicate_effect"]] <- rep(0, n_replicate)
    initial_values[["tau_replicate"]] <- 1

  }


  # ==========================================================================
  # 22. Temporal effect
  # ==========================================================================

  if (temporal_effect != "none") {
    if (is.null(temporal_data)) {
      stop(
        "temporal_data must be supplied when temporal_effect = '",
        temporal_effect,
        "'.",
        call. = FALSE
      )
    }
    if (length(temporal_data) != N) {
      stop(
        "temporal_data must have the same number of observations as ",
        "modelling_data.",
        call. = FALSE
      )
    }
    if (!is.numeric(temporal_data) || any(!is.finite(temporal_data))) {
      stop(
        "temporal_data must be numeric and contain only finite values.",
        call. = FALSE
      )
    }

    if (temporal_effect == "linear") {
      if (include_biomass) {
        initial_values[["beta_time_biomass"]] <- 0
      }
      if (include_diversity) {
        initial_values[["beta_time_diversity"]] <- 0
      }
      initial_values[["beta_time"]] <- 0

    }

    if (temporal_effect == "cubic_cyclic") {
      if (is.null(temporal_basis_penalized)) {
        stop(
          "temporal_basis_penalized must be supplied when ",
          "temporal_effect = 'cubic_cyclic'.",
          call. = FALSE
        )
      }
      if (!is.matrix(temporal_basis_penalized) ||
          !is.numeric(temporal_basis_penalized)) {
        stop(
          "temporal_basis_penalized must be a numeric matrix.",
          call. = FALSE
        )
      }
      if (nrow(temporal_basis_penalized) != N) {
        stop(
          "The number of rows in temporal_basis_penalized must ",
          "equal the number of observations in modelling_data.",
          call. = FALSE
        )
      }
      if (any(!is.finite(temporal_basis_penalized))) {
        stop(
          "temporal_basis_penalized contains non-finite values.",
          call. = FALSE
        )
      }

      n_temporal_basis <- ncol(temporal_basis_penalized)
      if (n_temporal_basis < 1) {
        stop(
          "temporal_basis_penalized contains no retained basis functions.",
          call. = FALSE
        )
      }

      if (isTRUE(temporal_shared)) {
        initial_values[["beta_time_shared"]] <- rep(
          0,
          n_temporal_basis
        )
        initial_values[["temporal_tau_shared"]] <- 1

        # When both responses are included, the model estimates a
        # diversity-specific loading; biomass defines the reference
        # scale with its loading fixed to one.
        if (include_biomass && include_diversity) {
          initial_values[["lambda_time_diversity"]] <- 1
        }
      } else {
        if (include_biomass) {
          initial_values[["beta_time_biomass"]] <- rep(
            0,
            n_temporal_basis
          )
          initial_values[["temporal_tau_biomass"]] <- 1
        }
        if (include_diversity) {
          initial_values[["beta_time_diversity"]] <- rep(
            0,
            n_temporal_basis
          )
          initial_values[["temporal_tau_diversity"]] <- 1
        }
      }
    }
  }


  # ==========================================================================
  # 23. Shared spatial effect
  # ==========================================================================

  if (spatial_effect == "coordinates") {
    initial_values[["beta_longitude_spatial"]] <- 0
    initial_values[["beta_latitude_spatial"]] <- 0
    initial_values[["beta_longitude_spatial_squared"]] <- 0
    initial_values[["beta_latitude_spatial_squared"]] <- 0
  }

  if (spatial_effect == "autocorrelation") {
    # Count unique projected locations using the same key convention as
    # the spatial-constants builder.
    spatial_points <- modelling_data |>
      dplyr::select(dplyr::all_of(c(longitude_col, latitude_col))) |>
      sf::st_as_sf(
        coords = c(longitude_col, latitude_col),
        crs = 4326,
        remove = FALSE
      ) |>
      sf::st_transform(crs = 25833)

    xy <- sf::st_coordinates(spatial_points)
    spatial_key <- paste(
      round(xy[, 1], 6),
      round(xy[, 2], 6),
      sep = "_"
    )
    n_spatial <- length(unique(spatial_key))

    if (n_spatial < 1) {
      stop("No unique spatial locations were identified.", call. = FALSE)
    }

    if (isTRUE(spatial_shared)) {
      initial_values[["spatial_effect"]] <- rep(0, n_spatial)
      initial_values[["rho_spatial"]] <- mean(spatial_range_prior)
      initial_values[["sigma_spatial"]] <- max(0.5, spatial_sigma_prior)

      # The spatiotemporal shared builder identifies biomass loading as 1
      # and estimates a diversity loading for the shared spatial field.
      if (model_structure == "spatiotemporal_shared" &&
          include_biomass && include_diversity) {
        initial_values[["lambda_spatial_diversity"]] <- 1
      }
    } else {
      if (include_biomass) {
        initial_values[["spatial_effect_biomass"]] <- rep(0, n_spatial)
        initial_values[["rho_spatial_biomass"]] <- mean(spatial_range_prior)
        initial_values[["sigma_spatial_biomass"]] <- max(0.5, spatial_sigma_prior)
        initial_values[["lambda_spatial_biomass"]] <- 0
      }
      if (include_diversity) {
        initial_values[["spatial_effect_diversity"]] <- rep(0, n_spatial)
        initial_values[["rho_spatial_diversity"]] <- mean(spatial_range_prior)
        initial_values[["sigma_spatial_diversity"]] <- max(0.5, spatial_sigma_prior)
        initial_values[["lambda_spatial_diversity"]] <- 0
      }
    }
  }


  # ==========================================================================
  # 24. Adapt shared-core values to the parameterisation of each builder
  # ==========================================================================

  # The baseline, separate temporal/spatial, and separate spatiotemporal
  # builders use response-specific fixed and ordinary random effects. The
  # shared temporal and spatial builders use shared fixed/random effects. The
  # spatiotemporal shared builder shares fixed effects but retains response-
  # specific year/locality/replicate effects in its current implementation.
  separate_fixed <- !share_fixed_effects
  separate_categorical <- separate_fixed
  separate_random <- !share_hyperparameters #separate_fixed || model_structure == "spatiotemporal_shared"

  if (separate_fixed && "beta" %in% names(initial_values)) {
    shared_beta <- initial_values[["beta"]]
    initial_values$beta <- NULL
    if (include_biomass) initial_values$beta_biomass <- shared_beta
    if (include_diversity) initial_values$beta_diversity <- shared_beta
  }

  if (separate_categorical) {
    # Categorical effects are shared in the spatiotemporal-shared builder,
    # but response-specific in the baseline and separate-effect builders.
    for (cov in categorical_covariates) {
      cov_clean <- make.names(cov)
      effect_name <- paste0(cov_clean, "_effect")
      tau_name <- paste0("tau_", cov_clean)
      if (effect_name %in% names(initial_values)) {
        z <- initial_values[[effect_name]]
        initial_values[[effect_name]] <- NULL
        if (include_biomass) initial_values[[paste0(cov_clean, "_effect_biomass")]] <- z
        if (include_diversity) initial_values[[paste0(cov_clean, "_effect_diversity")]] <- z
      }
      if (tau_name %in% names(initial_values)) {
        z <- initial_values[[tau_name]]
        initial_values[[tau_name]] <- NULL
        if (include_biomass) initial_values[[paste0(tau_name, "_biomass")]] <- z
        if (include_diversity) initial_values[[paste0(tau_name, "_diversity")]] <- z
      }
    }
  }

  if (separate_random) {
    for (component in c("year", "locality", "replicate")) {
      effect_name <- paste0(component, "_effect")
      tau_name <- paste0("tau_", component)
      if (effect_name %in% names(initial_values)) {
        z <- initial_values[[effect_name]]
        initial_values[[effect_name]] <- NULL
        if (include_biomass) initial_values[[paste0(effect_name, "_biomass")]] <- z
        if (include_diversity) initial_values[[paste0(effect_name, "_diversity")]] <- z
      }
      if (tau_name %in% names(initial_values)) {
        z <- initial_values[[tau_name]]
        initial_values[[tau_name]] <- NULL
        if (include_biomass) initial_values[[paste0(tau_name, "_biomass")]] <- z
        if (include_diversity) initial_values[[paste0(tau_name, "_diversity")]] <- z
      }
    }
  }

  # For spatiotemporal models the shared spatial field has a response loading;
  # the spatial-only shared builder uses the field directly with no loading.
  if (model_structure == "spatiotemporal_shared" &&
      "spatial_effect" %in% names(initial_values) &&
      include_biomass && include_diversity &&
      !("lambda_spatial_diversity" %in% names(initial_values))) {
    initial_values$lambda_spatial_diversity <- 1
  }

  # ==========================================================================
  # 25. Return initial values
  # ==========================================================================

  initial_values

}


#' Generate Initial Values for All Seven Insect NIMBLE Models
#'
#' Convenience wrapper that calls [create_insect_nimble_initial_values()] once
#' for each of the seven model structures and returns a named list of initial-
#' value lists. Arguments are forwarded to the single-model function, except
#' `model_structure`, `temporal_effect`, `temporal_shared`, `spatial_effect`, and
#' `spatial_shared`, which are selected internally for each model.
#'
#' For this all-model wrapper, temporal models use the cubic cyclic spline and
#' spatial models use spatial autocorrelation. Supply `temporal_data` and
#' `temporal_basis_penalized` for the temporal models, and valid longitude and
#' latitude columns in `modelling_data` for spatial models.
#'
#' @param ... Arguments forwarded to [create_insect_nimble_initial_values()],
#'   such as `modelling_data`, `modelCovs`, response columns, covariate types,
#'   random-effect flags, temporal basis data, and spatial priors.
#'
#' @return A named list with one element per model structure; each element is a
#'   named list of initial parameter values for that model.
#' @export
create_insect_nimble_initial_values_all_models <- function(...) {
  args <- list(...)

  # These arguments vary by builder and must be resolved per model.
  structure_specific <- c(
    "model_structure", "temporal_effect", "temporal_shared",
    "spatial_effect", "spatial_shared"
  )
  args[intersect(names(args), structure_specific)] <- NULL

  structures <- c(
    "baseline",
    "temporal_shared",
    "temporal_separate",
    "spatial_shared",
    "spatial_separate",
    "spatiotemporal_shared",
    "spatiotemporal_separate"
  )

  values <- lapply(structures, function(structure) {
    do.call(
      create_insect_nimble_initial_values,
      c(args, list(model_structure = structure))
    )
  })

  names(values) <- structures
  values
}
