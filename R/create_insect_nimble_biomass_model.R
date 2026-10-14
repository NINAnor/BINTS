#' Create a NIMBLE Model for Insect Biomass and Diversity
#'
#' Generate a Bayesian NIMBLE model specification for insect biomass,
#' species diversity, or both responses. The function acts as a dispatcher
#' that selects one of seven standalone model builders according to the
#' requested temporal and spatial structure.
#'
#' The available model structures are:
#' \describe{
#'   \item{Baseline}{No temporal or spatial process.}
#'   \item{Temporal shared}{A temporal process shared between biomass and
#'     diversity.}
#'   \item{Temporal separate}{Response-specific temporal processes.}
#'   \item{Spatial shared}{A spatial autocorrelation process shared between
#'     biomass and diversity.}
#'   \item{Spatial separate}{Response-specific spatial autocorrelation
#'     processes.}
#'   \item{Spatiotemporal shared}{Shared temporal and spatial processes.}
#'   \item{Spatiotemporal separate}{Separate temporal and spatial processes
#'     for the responses.}
#' }
#'
#' The selected model builder determines the implementation of the requested
#' effects. Spatiotemporal models combine temporal and spatial effects
#' additively; they do not include a location-specific temporal interaction.
#' The baseline model excludes both temporal and spatial processes but can
#' retain the optional grouping-level random effects.
#'
#' Continuous covariates are represented by fixed effects, while categorical
#' covariates can be represented by random effects according to the selected
#' model builder. Optional grouping-level random intercepts can be specified
#' for year, locality, and replicate.
#'
#' The function validates the model configuration, identifies the appropriate
#' standalone builder, and passes the relevant arguments to that builder.
#' The standalone model-builder functions must be defined in the calling
#' environment before this dispatcher is used.
#'
#' @param modelCovs Character vector specifying the names of covariates
#'   included in the model. The corresponding variables or design-matrix
#'   columns must be available in the model data supplied when fitting the
#'   generated NIMBLE model.
#'
#' @param response Character string specifying the response to model.
#'   Supported values are \code{"biomass"}, \code{"diversity"}, and
#'   \code{"both"}. Defaults to \code{"biomass"} through the first element
#'   of the argument choices.
#'
#' @param biomass_distribution Character string specifying the observation
#'   distribution for biomass. Supported values are \code{"normal"} and
#'   \code{"lognormal"}. Defaults to \code{"normal"}.
#'
#' @param diversity_distribution Character string specifying the observation
#'   distribution for diversity. Supported values are \code{"poisson"} and
#'   \code{"negative_binomial"}. Defaults to \code{"poisson"}.
#'
#' @param covariate_types Optional named character vector specifying whether
#'   each covariate is \code{"continuous"} or \code{"categorical"}. Names
#'   must correspond to elements of \code{modelCovs}. The default is
#'   \code{NULL}; interpretation of unspecified covariate types is delegated
#'   to the selected model builder.
#'
#' @param random_year Logical. If \code{TRUE}, include a year-level random
#'   intercept where supported by the selected model builder. Defaults to
#'   \code{FALSE}.
#'
#' @param random_locality Logical. If \code{TRUE}, include a locality-level
#'   random intercept where supported by the selected model builder. Defaults
#'   to \code{FALSE}.
#'
#' @param random_replicate Logical. If \code{TRUE}, include a replicate-level
#'   random intercept where supported by the selected model builder. Defaults
#'   to \code{FALSE}.
#'
#' @param temporal_effect Character string specifying the temporal structure.
#'   Supported values are \code{"none"}, \code{"linear"}, and
#'   \code{"cubic_cyclic"}. If \code{"none"} and \code{spatial_effect} is
#'   also \code{"none"}, the baseline model is selected. Otherwise, the
#'   temporal setting determines which temporal or spatiotemporal builder is
#'   selected. Defaults to \code{"none"}.
#'
#' @param temporal_basis Optional numeric matrix containing the temporal
#'   spline basis evaluated at the observations. Required by model builders
#'   that use a cyclic cubic spline. Its number of columns must be compatible
#'   with \code{temporal_penalty}.
#'
#' @param temporal_penalty Optional numeric square matrix representing the
#'   penalty associated with \code{temporal_basis}. Required by model builders
#'   that use a cyclic cubic spline. Its dimensions must match the number of
#'   columns in \code{temporal_basis}.
#'
#' @param temporal_penalty_tol Positive finite numeric scalar controlling
#'   the retention of positive eigenvalues when decomposing the temporal
#'   penalty matrix. Defaults to \code{1e-8}.
#'
#' @param temporal_shared Logical. Determines whether the temporal process is
#'   shared between biomass and diversity when the selected model structure
#'   includes a temporal effect. Shared temporal models currently require
#'   \code{temporal_effect = "cubic_cyclic"}. For spatiotemporal models, this
#'   setting must match \code{spatial_shared}. Defaults to \code{TRUE}.
#'
#' @param spatial_effect Character string specifying the spatial structure.
#'   Supported values are \code{"none"}, \code{"coordinates"}, and
#'   \code{"autocorrelation"}. If \code{"none"} and
#'   \code{temporal_effect} is also \code{"none"}, the baseline model is
#'   selected. Defaults to \code{"none"}.
#'
#' @param spatial_shared Logical. For spatial autocorrelation models, specifies
#'   whether the spatial process is shared between biomass and diversity.
#'   For spatiotemporal models, this setting must match
#'   \code{temporal_shared}. Separate coordinate-based spatial trends are not
#'   supported. Defaults to \code{TRUE}.
#'
#' @param spatial_range_prior Numeric vector of length two specifying the
#'   lower and upper bounds of the uniform prior for the spatial range
#'   parameter. Both values must be positive, and the upper bound must exceed
#'   the lower bound. Defaults to \code{c(1000, 500000)}.
#'
#' @param spatial_sigma_prior Positive finite numeric scalar specifying the
#'   scale used for the half-normal prior on the spatial standard deviation.
#'   Defaults to \code{1}.
#'
#' @param output_dir Character string specifying the directory in which the
#'   selected model builder writes the generated model script. The directory
#'   handling is delegated to that builder. Defaults to \code{temp_dir}.
#'
#' @param model_name Character string specifying the base name of the
#'   generated model script, without the \code{.R} extension. Defaults to
#'   \code{"biomass_diversity_model_nimble"}.
#'
#' @return Returns the object produced by the selected standalone model
#'   builder. This is expected to be a list containing the generated model
#'   code, output file path, response and covariate configuration, and
#'   temporal or spatial metadata. When the result is a list, the dispatcher
#'   adds the following elements:
#'   \describe{
#'     \item{selected_function}{Name of the standalone model-builder function
#'       selected by the dispatcher.}
#'     \item{effect_structure}{Selected structure: \code{"none"},
#'       \code{"temporal"}, \code{"spatial"}, or \code{"spatiotemporal"}.}
#'     \item{temporal_shared}{Requested temporal-sharing setting when a
#'       temporal effect is included; otherwise \code{NA}.}
#'     \item{spatial_shared}{Requested spatial-sharing setting when a
#'       spatial effect is included; otherwise \code{NA}.}
#'   }
#'
#' @details
#' The dispatcher selects a model builder according to the following rules:
#' \describe{
#'   \item{No temporal or spatial effect}{
#'     \code{create_insect_nimble_baseline_model()}.}
#'   \item{Temporal effect only}{
#'     \code{create_temporal_shared_model()} or
#'     \code{create_temporal_separate_model()}, according to
#'     \code{temporal_shared}.}
#'   \item{Spatial effect only}{
#'     \code{create_spatial_shared_model()} or
#'     \code{create_spatial_separate_model()}, according to
#'     \code{spatial_shared}.}
#'   \item{Both effects}{
#'     \code{create_spatiotemporal_shared_model()} or
#'     \code{create_spatiotemporal_separate_model()}. Both sharing settings
#'     must be \code{TRUE} or both must be \code{FALSE}.}
#' }
#'
#' The dispatcher does not itself construct or fit the NIMBLE model. It
#' delegates model generation to the selected function. The caller must
#' provide compatible data, constants, initial values, and MCMC settings
#' when compiling and fitting the generated model.
#'
#' The baseline excludes temporal and spatial processes but does not
#' automatically disable the optional year, locality, or replicate random
#' intercepts. Set the corresponding \code{random_*} arguments to
#' \code{FALSE} if those terms should also be excluded.
#'
#' Shared temporal models currently require a cyclic cubic spline.
#' Spatiotemporal models support shared temporal and shared spatial effects,
#' or separate temporal and separate spatial effects. Mixed sharing
#' configurations are not implemented by the available model builders.
#' Separate coordinate-based spatial trends are not supported.
#'
#' @examples
#' \dontrun{
#' # Baseline model: no temporal or spatial process.
#' baseline_spec <- create_insect_nimble_model(
#'   modelCovs = c("temperature", "habitat"),
#'   response = "both",
#'   biomass_distribution = "lognormal",
#'   diversity_distribution = "negative_binomial",
#'   temporal_effect = "none",
#'   spatial_effect = "none",
#'   random_year = FALSE,
#'   random_locality = TRUE,
#'   random_replicate = FALSE,
#'   output_dir = tempdir(),
#'   model_name = "insect_baseline"
#' )
#'
#' # Temporal model with a shared cyclic spline.
#' temporal_spec <- create_insect_nimble_model(
#'   modelCovs = c("temperature", "habitat"),
#'   response = "both",
#'   temporal_effect = "cubic_cyclic",
#'   temporal_basis = temporal_basis,
#'   temporal_penalty = temporal_penalty,
#'   temporal_shared = TRUE,
#'   spatial_effect = "none",
#'   output_dir = tempdir(),
#'   model_name = "insect_temporal_shared"
#' )
#'
#' # Spatial autocorrelation model with separate spatial processes.
#' spatial_spec <- create_insect_nimble_model(
#'   modelCovs = c("temperature", "habitat"),
#'   response = "both",
#'   temporal_effect = "none",
#'   spatial_effect = "autocorrelation",
#'   spatial_shared = FALSE,
#'   output_dir = tempdir(),
#'   model_name = "insect_spatial_separate"
#' )
#' }
#'
#' @seealso
#' \code{\link{create_insect_nimble_baseline_model}},
#' \code{\link{create_insect_temporal_shared_model}},
#' \code{\link{create_insect_temporal_separate_model}},
#' \code{\link{create_insect_spatial_shared_model}},
#' \code{\link{create_insect_spatial_separate_model}},
#' \code{\link{create_insect_spatiotemporal_shared_model}},
#' \code{\link{create_insect_spatiotemporal_separate_model}}
#'
#' @export

create_insect_nimble_model <- function(
    modelCovs,
    response = c("biomass", "diversity", "both"),
    biomass_distribution = c("normal", "lognormal"),
    diversity_distribution = c("poisson", "negative_binomial"),
    covariate_types = NULL,
    random_year = FALSE,
    random_locality = FALSE,
    random_replicate = FALSE,
    temporal_effect = c("none", "linear", "cubic_cyclic"),
    temporal_basis = NULL,
    temporal_penalty = NULL,
    temporal_penalty_tol = 1e-8,
    temporal_shared = TRUE,
    spatial_effect = c("none", "coordinates", "autocorrelation"),
    spatial_shared = TRUE,
    spatial_range_prior = c(1000, 500000),
    spatial_sigma_prior = 1,
    output_dir = temp_dir,
    model_name = "biomass_diversity_model_nimble"
) {
  # --------------------------------------------------------------------------
  # 1. Validate user-facing arguments
  # --------------------------------------------------------------------------
  response <- match.arg(response)
  biomass_distribution <- match.arg(biomass_distribution)
  diversity_distribution <- match.arg(diversity_distribution)
  temporal_effect <- match.arg(temporal_effect)
  spatial_effect <- match.arg(spatial_effect)

  if (!is.logical(temporal_shared) || length(temporal_shared) != 1L ||
      is.na(temporal_shared)) {
    stop("'temporal_shared' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(spatial_shared) || length(spatial_shared) != 1L ||
      is.na(spatial_shared)) {
    stop("'spatial_shared' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(temporal_penalty_tol) ||
      length(temporal_penalty_tol) != 1L ||
      !is.finite(temporal_penalty_tol) || temporal_penalty_tol <= 0) {
    stop("'temporal_penalty_tol' must be a positive finite number.",
         call. = FALSE)
  }

  if (!is.numeric(spatial_range_prior) || length(spatial_range_prior) != 2L ||
      any(!is.finite(spatial_range_prior)) || spatial_range_prior[1] <= 0 ||
      spatial_range_prior[2] <= spatial_range_prior[1]) {
    stop("'spatial_range_prior' must contain two positive values with upper > lower.",
         call. = FALSE)
  }

  if (!is.numeric(spatial_sigma_prior) || length(spatial_sigma_prior) != 1L ||
      !is.finite(spatial_sigma_prior) || spatial_sigma_prior <= 0) {
    stop("'spatial_sigma_prior' must be a positive finite number.",
         call. = FALSE)
  }

  has_temporal <- temporal_effect != "none"
  has_spatial <- spatial_effect != "none"

  # If both effects are absent, use the standalone baseline builder.
  # The baseline contains no temporal process and no spatial process.
  is_baseline <- !has_temporal && !has_spatial

  # The available spatio-temporal builders implement shared/shared and
  # separate/separate combinations. Mixed sharing (shared temporal + separate
  # spatial, or vice versa) requires an additional dedicated builder.
  if (has_temporal && has_spatial && temporal_shared != spatial_shared) {
    stop(
      "The available spatio-temporal functions support either both effects ",
      "shared (temporal_shared = TRUE and spatial_shared = TRUE) or both effects ",
      "separate (both FALSE). Mixed sharing is not implemented yet.",
      call. = FALSE
    )
  }

  # The current coordinate-based spatial trend uses the same coordinate
  # coefficients in both responses; separate coordinate coefficients are not
  # implemented by the six-function set.
  if (spatial_effect == "coordinates" && !spatial_shared) {
    stop(
      "Separate coordinate-trend parameters for each response are not implemented yet.",
      "Set spatial_shared = TRUE or use spatial_effect = 'autocorrelation' ",
      "with spatial_shared = FALSE.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------
  # 2. Select the relevant standalone function
  # --------------------------------------------------------------------------
  if (is_baseline) {
    selected_function_name <- "create_insect_baseline_model"
  } else if (has_temporal && has_spatial) {
    selected_function_name <- if (temporal_shared) {
      "create_insect_spatiotemporal_shared_model"
    } else {
      "create_insect_spatiotemporal_separate_model"
    }
  } else if (has_temporal) {
    selected_function_name <- if (temporal_shared) {
      "create_insect_temporal_shared_model"
    } else {
      "create_insect_temporal_separate_model"
    }
  } else {
    selected_function_name <- if (spatial_shared) {
      "create_insect_spatial_shared_model"
    } else {
      "create_insect_spatial_separate_model"
    }
  }

  if (!exists(selected_function_name, mode = "function", inherits = TRUE)) {
    stop(
      "Could not find '", selected_function_name, "()'. Source the standalone ",
      "model-builder function definitions before calling this dispatcher.",
      call. = FALSE
    )
  }

  selected_function <- get(selected_function_name, mode = "function", inherits = TRUE)

  # --------------------------------------------------------------------------
  # 3. Pass compatible arguments to the selected standalone function
  # --------------------------------------------------------------------------
  # Explicitly build the argument list to avoid forwarding dispatcher-only
  # arguments or accidentally passing unsupported arguments through ... .
  args <- list(
    modelCovs = modelCovs,
    response = response,
    biomass_distribution = biomass_distribution,
    diversity_distribution = diversity_distribution,
    covariate_types = covariate_types,
    random_year = random_year,
    random_locality = random_locality,
    random_replicate = random_replicate,
    temporal_effect = temporal_effect,
    temporal_basis = temporal_basis,
    temporal_penalty = temporal_penalty,
    temporal_penalty_tol = temporal_penalty_tol,
    temporal_shared = temporal_shared,
    spatial_effect = spatial_effect,
    spatial_shared = spatial_shared,
    spatial_range_prior = spatial_range_prior,
    spatial_sigma_prior = spatial_sigma_prior,
    output_dir = output_dir,
    model_name = model_name
  )

  # Each standalone function fixes its own effect structure. Remove arguments
  # that are irrelevant to that selected structure, so they cannot override
  # the configuration embedded in the standalone function.
  if (selected_function_name == "create_insect_baseline_model") {
    args$temporal_effect <- "none"
    args$temporal_basis <- NULL
    args$temporal_penalty <- NULL
    args$temporal_shared <- FALSE
    args$spatial_effect <- "none"
    args$spatial_shared <- TRUE
  } else if (selected_function_name %in% c(
    "create_insect_temporal_shared_model", "create_insect_temporal_separate_model"
  )) {
    args$spatial_effect <- "none"
    args$spatial_shared <- TRUE
  } else if (selected_function_name %in% c(
    "create_insect_spatial_shared_model", "create_insect_spatial_separate_model"
  )) {
    args$temporal_effect <- "none"
    args$temporal_basis <- NULL
    args$temporal_penalty <- NULL
    args$temporal_shared <- FALSE
  }

  result <- do.call(selected_function, args)

  # --------------------------------------------------------------------------
  # 4. Add dispatch metadata to the returned model specification
  # --------------------------------------------------------------------------
  if (is.list(result)) {
    result$selected_function <- selected_function_name
    result$effect_structure <- if (is_baseline) {
      "none"
    } else if (has_temporal && has_spatial) {
      "spatiotemporal"
    } else if (has_temporal) {
      "temporal"
    } else {
      "spatial"
    }
    result$temporal_shared <- if (has_temporal) temporal_shared else NA
    result$spatial_shared <- if (has_spatial) spatial_shared else NA
  }

  return(result)
}
