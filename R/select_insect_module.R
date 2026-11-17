
#' Select a Model Module for Insect Biomass and Diversity
#'
#' Selects and validates a model configuration specifying the temporal
#' effect, spatial effect, sharing structure, and random-effects structure
#' for an insect biomass and diversity model. The returned configuration
#' can be passed to [run_insect_nimble_model()] to construct and fit the
#' selected model.
#'
#' The function supports a baseline model, temporal-only models,
#' spatial-only models, and combined spatiotemporal models. For spatial
#' and spatiotemporal modules, the locality random effect is disabled to
#' avoid including an additional locality-level effect alongside the
#' selected spatial structure.
#'
#' Temporal and spatial effects can be shared between the modelled
#' responses where supported. Shared linear temporal effects and separate
#' coordinate-based spatial effects are not currently supported.
#'
#' @param module Character string specifying the model module to select.
#'   Supported values are:
#'   \describe{
#'     \item{`NULL`}{Selects the baseline model without an explicit temporal
#'       or spatial effect. Year, locality, and replicate random effects
#'       remain enabled. The character string `"NULL"` is also accepted.}
#'     \item{`"temporal_linear"`}{Includes a linear temporal effect without
#'       an explicit spatial effect. The temporal effect must be separate.}
#'     \item{`"temporal_cyclic"`}{Includes a cyclic cubic spline temporal
#'       effect without an explicit spatial effect. The temporal effect
#'       can be shared or separate.}
#'     \item{`"spatial_coordinates"`}{Selects the shared coordinate-based
#'       spatial effect without an explicit temporal effect.}
#'     \item{`"spatial_autocorrelation"`}{Selects a spatial autocorrelation
#'       effect without an explicit temporal effect. The spatial effect
#'       can be shared or separate.}
#'     \item{`"spatio_temporal_linear_coordinates"`}{Combines a linear
#'       temporal effect with a coordinate-based spatial effect. This
#'       combination is incompatible with the current sharing constraints
#'       because linear temporal effects must be separate and coordinate-
#'       based spatial effects must be shared.}
#'     \item{`"spatio_temporal_linear_autocorrelation"`}{Combines a linear
#'       temporal effect with spatial autocorrelation. Both effects must
#'       be separate under the current implementation.}
#'     \item{`"spatio_temporal_cyclic_coordinates"`}{Combines a cyclic cubic
#'       spline temporal effect with a coordinate-based spatial effect.
#'       Both effects must be shared under the current implementation.}
#'     \item{`"spatio_temporal_cyclic_autocorrelation"`}{Combines a cyclic
#'       cubic spline temporal effect with spatial autocorrelation. Both
#'       effects can be shared or separate, but their sharing settings
#'       must agree.}
#'   }
#'   Defaults to `NULL`.
#'
#' @param temporal_basis_k Integer specifying the requested basis dimension
#'   for the cyclic cubic spline temporal effect. Must be a finite integer
#'   greater than or equal to four. Defaults to `10`. This argument is
#'   retained in the returned configuration for all modules, but is used
#'   to construct the temporal spline only when `temporal_effect` is
#'   `"cubic_cyclic"`.
#'
#' @param temporal_shared Optional logical value specifying whether the
#'   temporal effect is shared between the modelled responses. Must be a
#'   single non-missing logical value when supplied. If `NULL`, defaults
#'   to `TRUE` for a cyclic cubic temporal effect and `FALSE` otherwise.
#'   Shared temporal effects currently require
#'   `temporal_effect = "cubic_cyclic"`. Must be `FALSE` when
#'   `temporal_effect = "none"`.
#'
#' @param spatial_shared Optional logical value specifying whether the
#'   spatial effect is shared between the modelled responses. Must be a
#'   single non-missing logical value when supplied. If `NULL`, defaults
#'   to `FALSE` when no spatial effect is selected, `TRUE` for a
#'   spatial-only model, and the value of `temporal_shared` for a
#'   spatiotemporal model. Must be `FALSE` when `spatial_effect = "none"`.
#'   Separate coordinate-based spatial effects are not supported.
#'   Spatiotemporal models require `temporal_shared` and `spatial_shared`
#'   to have the same value.
#'
#' @details
#' The function returns a named list containing the following
#' configuration elements:
#' \describe{
#'   \item{`random_year`}{Logical indicating whether year-level random
#'     effects are enabled. Set to `TRUE` for all modules.}
#'   \item{`random_locality`}{Logical indicating whether locality-level
#'     random effects are enabled. Set to `TRUE` for the baseline and
#'     temporal-only modules, and `FALSE` for spatial and spatiotemporal
#'     modules.}
#'   \item{`random_replicate`}{Logical indicating whether replicate-level
#'     random effects are enabled. Set to `TRUE` for all modules.}
#'   \item{`temporal_effect`}{Character string specifying `"none"`,
#'     `"linear"`, or `"cubic_cyclic"`.}
#'   \item{`temporal_shared`}{Logical indicating whether the temporal
#'     effect is shared between modelled responses.}
#'   \item{`temporal_basis_k`}{Integer specifying the requested cyclic
#'     spline basis dimension.}
#'   \item{`spatial_effect`}{Character string specifying `"none"`,
#'     `"coordinates"`, or `"autocorrelation"`.}
#'   \item{`spatial_shared`}{Logical indicating whether the spatial
#'     effect is shared between modelled responses.}
#' }
#'
#' The selected configuration determines model structure but does not
#' supply the observation data or spatial and temporal variables. Those
#' inputs must be supplied separately to [run_insect_nimble_model()].
#' Temporal models require the columns specified by `temporal_col`;
#' spatial models require the coordinate columns specified by
#' `longitude_col` and `latitude_col`. Spatial autocorrelation models
#' additionally require appropriate spatial range and spatial
#' standard-deviation priors when the model is fitted.
#'
#' The function prints a message describing the selected module and its
#' final configuration. The returned configuration is invisible.
#'
#' @return A named list containing the model configuration. The list can
#'   be supplied as arguments to [run_insect_nimble_model()] using
#'   `do.call()`.
#'
#' @seealso
#' [run_insect_nimble_model()],
#' [create_insect_nimble_model()]
#'
#' @examples
#' # Baseline model
#' baseline_config <- select_insect_module()
#'
#' # Linear temporal effect: separate by default
#' linear_config <- select_insect_module("temporal_linear")
#'
#' # Shared cyclic cubic temporal effect
#' cyclic_config <- select_insect_module(
#'   "temporal_cyclic",
#'   temporal_basis_k = 10
#' )
#'
#' # Separate cyclic cubic temporal effects
#' cyclic_separate_config <- select_insect_module(
#'   "temporal_cyclic",
#'   temporal_shared = FALSE
#' )
#'
#' # Shared spatial autocorrelation
#' spatial_shared_config <- select_insect_module(
#'   "spatial_autocorrelation"
#' )
#'
#' # Separate spatial autocorrelation
#' spatial_separate_config <- select_insect_module(
#'   "spatial_autocorrelation",
#'   spatial_shared = FALSE
#' )
#'
#' # Shared cyclic spatiotemporal model
#' cyclic_spatiotemporal_config <- select_insect_module(
#'   "spatio_temporal_cyclic_autocorrelation"
#' )
#'
#' # Separate cyclic spatiotemporal model
#' cyclic_spatiotemporal_separate_config <- select_insect_module(
#'   "spatio_temporal_cyclic_autocorrelation",
#'   temporal_shared = FALSE,
#'   spatial_shared = FALSE
#' )
#'
#' # Linear spatiotemporal autocorrelation model:
#' # both effects are separate under the current constraints
#' linear_spatiotemporal_config <- select_insect_module(
#'   "spatio_temporal_linear_autocorrelation"
#' )
#'
#' # Pass a selected configuration to the model-fitting function
#' # fit <- do.call(
#' #   run_insect_nimble_model,
#' #   c(
#' #     list(
#' #       modelling_data = modelling_data,
#' #       modelCovs = modelCovs,
#' #       covariate_types = covariate_types
#' #     ),
#' #     cyclic_spatiotemporal_config
#' #   )
#' # )
#'
#' @export



select_insect_module <- function(
    module = NULL,
    temporal_basis_k = 10,
    temporal_shared = NULL,
    spatial_shared = NULL
) {

  # --------------------------------------------------------
  # 1. Handle NULL input
  # --------------------------------------------------------

  if (is.character(module) &&
      length(module) == 1L &&
      !is.na(module) &&
      module == "NULL") {
    module <- NULL
  }

  # --------------------------------------------------------
  # 2. Validate module
  # --------------------------------------------------------

  valid_modules <- c(
    "baseline",
    "temporal_linear",
    "temporal_cyclic",
    "spatial_coordinates",
    "spatial_autocorrelation",
    "spatio_temporal_linear_coordinates",
    "spatio_temporal_linear_autocorrelation",
    "spatio_temporal_cyclic_coordinates",
    "spatio_temporal_cyclic_autocorrelation"
  )

  if (!is.null(module) &&
      (!is.character(module) ||
       length(module) != 1L ||
       is.na(module) ||
       !module %in% valid_modules)) {

    stop(
      "'module' must be NULL or one of: ",
      paste(valid_modules, collapse = ", "),
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # 3. Validate temporal basis dimension
  # --------------------------------------------------------

  if (!is.numeric(temporal_basis_k) ||
      length(temporal_basis_k) != 1L ||
      is.na(temporal_basis_k) ||
      !is.finite(temporal_basis_k) ||
      temporal_basis_k < 4 ||
      temporal_basis_k != as.integer(temporal_basis_k)) {

    stop(
      "'temporal_basis_k' must be a single integer >= 4.",
      call. = FALSE
    )
  }

  temporal_basis_k <- as.integer(temporal_basis_k)

  # --------------------------------------------------------
  # 4. Validate sharing arguments
  # --------------------------------------------------------

  validate_shared_flag <- function(x, arg_name) {

    if (!is.null(x) &&
        (!is.logical(x) ||
         length(x) != 1L ||
         is.na(x))) {

      stop(
        "'", arg_name,
        "' must be NULL, TRUE, or FALSE.",
        call. = FALSE
      )
    }
  }

  validate_shared_flag(
    temporal_shared,
    "temporal_shared"
  )

  validate_shared_flag(
    spatial_shared,
    "spatial_shared"
  )

  # --------------------------------------------------------
  # 5. Initialise configuration
  # --------------------------------------------------------

  config <- list(
    random_year = TRUE,
    random_locality = TRUE,
    random_replicate = TRUE,
    temporal_effect = "none",
    temporal_shared = FALSE,
    temporal_basis_k = temporal_basis_k,
    spatial_effect = "none",
    spatial_shared = FALSE
  )

  # --------------------------------------------------------
  # 6. Define model components
  # --------------------------------------------------------

  if (is.null(module)) {

    config_name <- "null"

  } else if (module == "baseline"){

    config_name <- "baseline"
    config$random_locality <- FALSE

  } else if (module == "temporal_linear") {

    config$temporal_effect <- "linear"
    config$random_locality <- FALSE
    config_name <- module

  } else if (module == "temporal_cyclic") {

    config$temporal_effect <- "cubic_cyclic"
    config$random_replicate <- FALSE
    config_name <- module

  } else if (module == "spatial_coordinates") {

    config$spatial_effect <- "coordinates"
    config$random_locality <- FALSE
      config$random_locality <- FALSE
    config_name <- module

  } else if (module == "spatial_autocorrelation") {

    config$spatial_effect <- "autocorrelation"
    config$random_locality <- FALSE
    config_name <- module

  } else if (module == "spatio_temporal_linear_coordinates") {

    config$temporal_effect <- "linear"
    config$spatial_effect <- "coordinates"
    config$random_locality <- FALSE
    config_name <- module

  } else if (module == "spatio_temporal_linear_autocorrelation") {

    config$temporal_effect <- "linear"
    config$spatial_effect <- "autocorrelation"
    config$random_locality <- FALSE
    config_name <- module

  } else if (module == "spatio_temporal_cyclic_coordinates") {

    config$temporal_effect <- "cubic_cyclic"
    config$spatial_effect <- "coordinates"
    config$random_locality <- FALSE
    config_name <- module

  } else if (module == "spatio_temporal_cyclic_autocorrelation") {

    config$temporal_effect <- "cubic_cyclic"
    config$spatial_effect <- "autocorrelation"
    config$random_locality <- FALSE
    config_name <- module
  }

  # --------------------------------------------------------
  # 7. Set default sharing configuration
  # --------------------------------------------------------

  # Linear temporal effects must be separate.
  # Cyclic cubic temporal effects can be shared or separate.
  if (is.null(temporal_shared)) {

    config$temporal_shared <-
      identical(config$temporal_effect, "cubic_cyclic")

  } else {

    config$temporal_shared <- temporal_shared
  }

  # Default spatial sharing:
  # - No spatial effect: FALSE
  # - Spatial-only model: TRUE
  # - Spatiotemporal model: follow temporal sharing
  if (is.null(spatial_shared)) {

    if (identical(config$spatial_effect, "none")) {

      config$spatial_shared <- FALSE

    } else if (!identical(config$temporal_effect, "none")) {

      config$spatial_shared <- config$temporal_shared

    } else {

      config$spatial_shared <- TRUE
    }

  } else {

    config$spatial_shared <- spatial_shared
  }

  # --------------------------------------------------------
  # 8. Validate model compatibility
  # --------------------------------------------------------

  # Shared linear temporal effects are not implemented.
  # if (identical(config$temporal_effect, "linear") &&
  #     config$temporal_shared) {
  #
  #   stop(
  #     "Shared temporal effects currently require ",
  #     "temporal_effect = 'cubic_cyclic'. ",
  #     "For a linear temporal effect, set temporal_shared = FALSE.",
  #     call. = FALSE
  #   )
  # }

  # Sharing is not meaningful when an effect is inactive.
  if (identical(config$temporal_effect, "none") &&
      config$temporal_shared) {

    stop(
      "temporal_shared must be FALSE when temporal_effect = 'none'.",
      call. = FALSE
    )
  }

  if (identical(config$spatial_effect, "none") &&
      config$spatial_shared) {

    stop(
      "spatial_shared must be FALSE when spatial_effect = 'none'.",
      call. = FALSE
    )
  }

  # Separate coordinate-based spatial trends are not implemented.
  # if (identical(config$spatial_effect, "coordinates") &&
  #     !config$spatial_shared) {
  #
  #   stop(
  #     "Separate coordinate-based spatial effects are not supported. ",
  #     "Use spatial_shared = TRUE or choose ",
  #     "spatial_effect = 'autocorrelation'.",
  #     call. = FALSE
  #   )
  # }

  # Spatiotemporal models require matching sharing structures.
  has_temporal <- !identical(
    config$temporal_effect,
    "none"
  )

  has_spatial <- !identical(
    config$spatial_effect,
    "none"
  )
  #
  # if (has_temporal &&
  #     has_spatial &&
  #     !identical(
  #       config$temporal_shared,
  #       config$spatial_shared
  #     )) {
  #
  #   stop(
  #     "Spatiotemporal models require temporal_shared and ",
  #     "spatial_shared to have the same value.",
  #     call. = FALSE
  #   )
  # }

  # --------------------------------------------------------
  # 9. Report selected configuration
  # --------------------------------------------------------

  message(
    "\nSelected module: ", config_name,
    "\n  Random year effect: ", config$random_year,
    "\n  Random locality effect: ", config$random_locality,
    "\n  Random replicate effect: ", config$random_replicate,
    "\n  Temporal effect: ", config$temporal_effect,
    "\n  Temporal shared: ", config$temporal_shared,
    "\n  Temporal basis dimension: ", config$temporal_basis_k,
    "\n  Spatial effect: ", config$spatial_effect,
    "\n  Spatial shared: ", config$spatial_shared,
    "\n"
  )

  invisible(config)
}
