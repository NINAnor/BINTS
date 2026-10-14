
#' Select a Model Module for Insect Biomass and Diversity
#'
#' Selects a model configuration specifying the temporal effect, spatial
#' effect, and random-effects structure for an insect biomass and diversity
#' model. The returned configuration can be passed to
#' [run_insect_nimble_model()] to construct and fit the selected model.
#'
#' The function supports a baseline model, temporal-only models, spatial-only
#' models, and combined spatiotemporal models. For spatial and spatiotemporal
#' modules, the locality random effect is disabled to avoid including an
#' additional locality-level effect alongside the selected spatial structure.
#'
#' @param module Character string specifying the model module to select.
#'   Supported values are:
#'   \describe{
#'     \item{`NULL`}{Selects the baseline model without an explicit temporal
#'       or spatial effect. Year, locality, and replicate random effects
#'       remain enabled. The character string `"NULL"` is also accepted.}
#'     \item{`"temporal_linear"`}{Includes a linear temporal effect without
#'       an explicit spatial effect.}
#'     \item{`"temporal_cyclic"`}{Includes a cyclic cubic spline temporal
#'       effect without an explicit spatial effect.}
#'     \item{`"spatial_coordinates"`}{Selects the coordinate-based spatial
#'       effect without an explicit temporal effect.}
#'     \item{`"spatial_autocorrelation"`}{Selects the spatial autocorrelation
#'       effect without an explicit temporal effect.}
#'     \item{`"spatio_temporal_linear_coordinates"`}{Combines a linear
#'       temporal effect with a coordinate-based spatial effect.}
#'     \item{`"spatio_temporal_linear_autocorrelation"`}{Combines a linear
#'       temporal effect with spatial autocorrelation.}
#'     \item{`"spatio_temporal_cyclic_coordinates"`}{Combines a cyclic cubic
#'       spline temporal effect with a coordinate-based spatial effect.}
#'     \item{`"spatio_temporal_cyclic_autocorrelation"`}{Combines a cyclic
#'       cubic spline temporal effect with spatial autocorrelation.}
#'   }
#'   Defaults to `NULL`.
#' @param temporal_basis_k Integer specifying the requested basis dimension
#'   for the cyclic cubic spline temporal effect. Must be a finite integer
#'   greater than or equal to four. Defaults to `10`. This argument is
#'   retained in the returned configuration for all modules, but is used
#'   to construct the temporal spline only when `temporal_effect` is
#'   `"cubic_cyclic"`.
#'
#' @details
#' The function returns a named list containing the following configuration
#' elements:
#' \describe{
#'   \item{`random_year`}{Logical indicating whether year-level random effects
#'     are enabled. Set to `TRUE` for all modules.}
#'   \item{`random_locality`}{Logical indicating whether locality-level random
#'     effects are enabled. Set to `TRUE` for the baseline and temporal-only
#'     modules, and `FALSE` for spatial and spatiotemporal modules.}
#'   \item{`random_replicate`}{Logical indicating whether replicate-level
#'     random effects are enabled. Set to `TRUE` for all modules.}
#'   \item{`temporal_effect`}{Character string passed to the model-fitting
#'     function. Takes one of `"none"`, `"linear"`, or `"cubic_cyclic"`.}
#'   \item{`temporal_basis_k`}{Integer specifying the requested cyclic spline
#'     basis dimension.}
#'   \item{`spatial_effect`}{Character string passed to the model-fitting
#'     function. Takes one of `"none"`, `"coordinates"`, or
#'     `"autocorrelation"`.}
#'   \item{`spatial_shared`}{Logical indicating whether the spatial structure
#'     is shared between modelled responses, where supported. Set to `TRUE`
#'     for all modules.}
#' }
#'
#' The selected configuration determines model structure but does not supply
#' the observation data or spatial and temporal variables. Those inputs must
#' be supplied separately to [run_insect_nimble_model()]. Temporal models
#' require the columns specified by `temporal_col`; spatial models require
#' the coordinate columns specified by `longitude_col` and `latitude_col`.
#' Spatial autocorrelation models may additionally use `spatial_range_prior`
#' and `spatial_sigma_prior`.
#'
#' Each call prints a message describing the selected module and its
#' configuration. The function also prints the final configuration before
#' returning it.
#'
#' @return A named list containing the model configuration. The list can be
#'   supplied as arguments to [run_insect_nimble_model()] using `do.call()`.
#'
#' @seealso
#' [run_insect_nimble_model()],
#' [create_insect_nimble_model()]
#'
#' @examples
#' # Baseline model
#' baseline_config <- select_insect_module()
#'
#' # Temporal-only models
#' linear_config <- select_insect_module("temporal_linear")
#' cyclic_config <- select_insect_module(
#'   "temporal_cyclic",
#'   temporal_basis_k = 10
#' )
#'
#' # Spatial-only models
#' coordinates_config <- select_insect_module("spatial_coordinates")
#' autocorrelation_config <- select_insect_module("spatial_autocorrelation")
#'
#' # Spatiotemporal models
#' linear_coordinates_config <- select_insect_module(
#'   "spatio_temporal_linear_coordinates"
#' )
#' linear_autocorrelation_config <- select_insect_module(
#'   "spatio_temporal_linear_autocorrelation"
#' )
#' cyclic_coordinates_config <- select_insect_module(
#'   "spatio_temporal_cyclic_coordinates",
#'   temporal_basis_k = 10
#' )
#' cyclic_autocorrelation_config <- select_insect_module(
#'   "spatio_temporal_cyclic_autocorrelation",
#'   temporal_basis_k = 10
#' )
#'
#' # Pass a selected configuration to the model-fitting function
#' # fit <- do.call(
#' #   run_insect_nimble_model,
#' #   c(
#' #     list(
#' #       modelling_data = modelling_data,
#' #       modelCovs = modelCovs
#' #     ),
#' #     cyclic_autocorrelation_config
#' #   )
#' # )
#'
#' @export

select_insect_module <- function(module = NULL, temporal_basis_k = 10) {

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
  # 4. Default configuration
  # --------------------------------------------------------

  config <- list(
    random_year = TRUE,
    random_locality = TRUE,
    random_replicate = TRUE,
    temporal_effect = "none",
    temporal_basis_k = temporal_basis_k,
    spatial_effect = "none",
    spatial_shared = TRUE
  )

  # --------------------------------------------------------
  # 5. Select model components
  # --------------------------------------------------------

  if (is.null(module)) {

    message(
      "\nModule: baseline model\n",
      "  Temporal effect: none\n",
      "  Spatial effect: none\n",
      "  Random effects: year, locality, and replicate\n"
    )

  } else if (identical(module, "temporal_linear")) {

    config$temporal_effect <- "linear"

    message(
      "\nModule: temporal_linear\n",
      "  Temporal effect: linear trend\n",
      "  Spatial effect: none\n",
      "  Input: temporal_data, derived from temporal_col\n",
      "  temporal_basis_k: retained but not used by the linear effect\n",
      "  Random effects: year, locality, and replicate\n"
    )

  } else if (identical(module, "temporal_cyclic")) {

    config$temporal_effect <- "cubic_cyclic"

    message(
      "\nModule: temporal_cyclic\n",
      "  Temporal effect: cyclic cubic spline\n",
      "  Spatial effect: none\n",
      "  Inputs: temporal_data and temporal_basis_k\n",
      "  temporal_data: derived from temporal_col\n",
      "  temporal_basis_k: spline basis dimension (k = ",
      temporal_basis_k, ")\n",
      "  Random effects: year, locality, and replicate\n"
    )

  } else if (identical(module, "spatial_coordinates")) {

    config$spatial_effect <- "coordinates"
    config$random_locality <- FALSE

    message(
      "\nModule: spatial_coordinates\n",
      "  Temporal effect: none\n",
      "  Spatial effect: coordinate-based\n",
      "  Inputs: longitude_col and latitude_col\n",
      "  Required data: numeric spatial coordinates\n",
      "  random_locality: FALSE to avoid an additional locality effect\n",
      "  spatial_shared: ", config$spatial_shared, "\n"
    )

  } else if (identical(module, "spatial_autocorrelation")) {

    config$spatial_effect <- "autocorrelation"
    config$random_locality <- FALSE

    message(
      "\nModule: spatial_autocorrelation\n",
      "  Temporal effect: none\n",
      "  Spatial effect: spatial autocorrelation\n",
      "  Inputs: longitude_col and latitude_col\n",
      "  Required data: numeric spatial coordinates\n",
      "  Additional inputs: spatial_range_prior and spatial_sigma_prior\n",
      "  random_locality: FALSE to avoid an additional locality effect\n",
      "  spatial_shared: ", config$spatial_shared, "\n"
    )

  } else if (identical(module, "spatio_temporal_linear_coordinates")) {

    config$temporal_effect <- "linear"
    config$spatial_effect <- "coordinates"
    config$random_locality <- FALSE

    message(
      "\nModule: spatio_temporal_linear_coordinates\n",
      "  Temporal effect: linear trend\n",
      "  Spatial effect: coordinate-based\n",
      "  Temporal input: temporal_col\n",
      "  Spatial inputs: longitude_col and latitude_col\n",
      "  temporal_basis_k: retained but not used by the linear effect\n",
      "  random_locality: FALSE\n"
    )

  } else if (identical(module, "spatio_temporal_linear_autocorrelation")) {

    config$temporal_effect <- "linear"
    config$spatial_effect <- "autocorrelation"
    config$random_locality <- FALSE

    message(
      "\nModule: spatio_temporal_linear_autocorrelation\n",
      "  Temporal effect: linear trend\n",
      "  Spatial effect: spatial autocorrelation\n",
      "  Temporal input: temporal_col\n",
      "  Spatial inputs: longitude_col and latitude_col\n",
      "  Additional inputs: spatial_range_prior and spatial_sigma_prior\n",
      "  temporal_basis_k: retained but not used by the linear effect\n",
      "  random_locality: FALSE\n"
    )

  } else if (identical(module, "spatio_temporal_cyclic_coordinates")) {

    config$temporal_effect <- "cubic_cyclic"
    config$spatial_effect <- "coordinates"
    config$random_locality <- FALSE

    message(
      "\nModule: spatio_temporal_cyclic_coordinates\n",
      "  Temporal effect: cyclic cubic spline\n",
      "  Spatial effect: coordinate-based\n",
      "  Temporal inputs: temporal_col and temporal_basis_k\n",
      "  Spatial inputs: longitude_col and latitude_col\n",
      "  temporal_basis_k: ", temporal_basis_k, "\n",
      "  random_locality: FALSE\n"
    )

  } else if (identical(module, "spatio_temporal_cyclic_autocorrelation")) {

    config$temporal_effect <- "cubic_cyclic"
    config$spatial_effect <- "autocorrelation"
    config$random_locality <- FALSE

    message(
      "\nModule: spatio_temporal_cyclic_autocorrelation\n",
      "  Temporal effect: cyclic cubic spline\n",
      "  Spatial effect: spatial autocorrelation\n",
      "  Temporal inputs: temporal_col and temporal_basis_k\n",
      "  Spatial inputs: longitude_col and latitude_col\n",
      "  Additional inputs: spatial_range_prior and spatial_sigma_prior\n",
      "  temporal_basis_k: ", temporal_basis_k, "\n",
      "  random_locality: FALSE\n"
    )
  }

  # --------------------------------------------------------
  # 6. Report final configuration
  # --------------------------------------------------------

  message(
    "\nSelected configuration:",
    "\n  random_year = ", config$random_year,
    "\n  random_locality = ", config$random_locality,
    "\n  random_replicate = ", config$random_replicate,
    "\n  temporal_effect = '", config$temporal_effect, "'",
    "\n  temporal_basis_k = ", config$temporal_basis_k,
    "\n  spatial_effect = '", config$spatial_effect, "'",
    "\n  spatial_shared = ", config$spatial_shared,
    "\n"
  )

  return(config)
}
