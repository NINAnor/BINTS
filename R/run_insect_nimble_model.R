#' Fit a NIMBLE Model for Insect Biomass and Diversity
#'
#' Prepare, build, compile, and optionally fit a Bayesian NIMBLE model for
#' insect biomass, diversity, or both responses. The function coordinates
#' model-code generation, data preparation, constant construction, initial
#' value generation, and NIMBLE model compilation.
#'
#' Optionally, the function configures and compiles an MCMC sampler, runs the
#' requested number of iterations, and calculates the Widely Applicable
#' Information Criterion (WAIC) when enabled.
#'
#' The model can include continuous covariates, categorical covariates treated
#' as random effects, year and locality random effects, replicate-level random
#' effects, temporal effects, and spatial effects.
#'
#' @param modelling_data A data frame containing the response variables,
#'   covariates, temporal columns, and any grouping or spatial coordinate
#'   columns required by the selected model.
#' @param modelCovs Character vector of covariate column names to include in
#'   the model. All specified columns must exist in `modelling_data`.
#' @param covariate_types Optional named character vector identifying each
#'   covariate as `"continuous"` or `"categorical"`. Every covariate in
#'   `modelCovs` must have a corresponding type. If `NULL`, all covariates are
#'   treated as continuous.
#' @param share_fixed_effects Logical; passed to the model-code and
#'   initial-value helpers. Defaults to `FALSE`.
#' @param share_hyperparameters Logical; passed to the model-code and
#'   initial-value helpers. Defaults to `FALSE`.
#' @param response Character string specifying the response to model:
#'   `"biomass"`, `"diversity"`, or `"both"`. Defaults to `"both"`.
#' @param biomass_col Name of the column containing the biomass response.
#'   Defaults to `"sum_wet_weight"`.
#' @param diversity_col Name of the column containing the diversity response.
#'   Defaults to `"shannon_div"`.
#' @param biomass_distribution Biomass likelihood: `"normal"` or
#'   `"lognormal"`. Defaults to `"normal"`.
#' @param diversity_distribution Diversity likelihood: `"normal"`, `"poisson"`,
#'   or `"negative_binomial"`. Defaults to `"normal"`.
#' @param replicate_col Name of the column containing replicate identifiers.
#'   Defaults to `"sampling_name"`. When supplied, the function constructs
#'   replicate indices using replicate numbers extracted from the suffix
#'   following the last hyphen, grouped by locality and year.
#' @param temporal_col Character vector of column names used to construct the
#'   temporal variable. Defaults to `c("start_date_julian", "end_date_julian")`.
#'   When a temporal effect is included, the function calculates the row mean
#'   of these columns.
#' @param random_year Logical; if `TRUE`, include year-level random effects.
#'   Defaults to `TRUE`.
#' @param random_locality Logical; if `TRUE`, include locality-level random
#'   effects. Defaults to `TRUE`.
#' @param random_replicate Logical; if `TRUE`, include replicate-level random
#'   effects. Defaults to `TRUE`.
#' @param temporal_effect Temporal effect specification: `"none"`, `"linear"`,
#'   or `"cubic_cyclic"`. Defaults to `"cubic_cyclic"`.
#' @param temporal_basis_k Integer specifying the requested basis dimension
#'   for the cubic cyclic spline. Must be an integer greater than or equal to
#'   four. Defaults to `10`.
#' @param temporal_shared Logical; for temporal effects, determines whether
#'   the temporal process is shared between biomass and diversity. Defaults to
#'   `TRUE`.
#' @param spatial_effect Spatial effect specification: `"none"`,
#'   `"coordinates"`, or `"autocorrelation"`. Defaults to `"none"`.
#' @param spatial_shared Logical; for spatial effects, determines whether the
#'   spatial structure is shared between the included responses. Defaults to
#'   `TRUE`.
#' @param longitude_col Name of the column containing longitude coordinates.
#'   Defaults to `"longitude"`.
#' @param latitude_col Name of the column containing latitude coordinates.
#'   Defaults to `"latitude"`.
#' @param spatial_range_prior Numeric vector of length two defining the lower
#'   and upper bounds of the spatial range prior, in metres. Defaults to
#'   `c(1000, 500000)`.
#' @param spatial_sigma_prior Positive finite numeric value specifying the
#'   spatial standard-deviation prior parameter. Defaults to `1`.
#' @param output_dir Directory in which the generated NIMBLE model source file
#'   is written. Defaults to `temp_dir`.
#' @param seed Integer seed(s) controlling reproducibility. Either a single
#'   value or a vector of length `nchains`. If a single value is supplied, the
#'   per-chain seeds are `seed, seed + 1, ..., seed + nchains - 1`. Each chain
#'   seed is used both to generate that chain's initial values and as the
#'   R random-number seed at the start of that chain in
#'   [nimble::runMCMC()]. Defaults to `123`.
#' @param run_mcmc Logical; if `TRUE`, configure, compile, and run MCMC after
#'   compiling the NIMBLE model. Defaults to `TRUE`.
#' @param run_waic Logical; if `TRUE`, enable WAIC calculation. Defaults to
#'   `TRUE`.
#' @param niter Total number of MCMC iterations requested. Defaults to `20000`.
#' @param nburnin Number of initial MCMC iterations to discard. Defaults to
#'   `5000`.
#' @param thin Thinning interval for the MCMC output. Defaults to `10`.
#' @param nchains Number of MCMC chains to run. Defaults to `2`.
#'
#' @details
#' The function carries out the following steps:
#' \enumerate{
#'   \item Validates the response specification, response distributions,
#'   covariate types, modelling data, spatial settings, seeds and chains.
#'   \item Extracts the selected response variables and constructs replicate
#'   indices when a replicate column is supplied.
#'   \item Constructs the temporal variable from the row means of the columns
#'   specified by `temporal_col`, when a temporal effect is requested.
#'   \item Standardises the temporal variable for a linear temporal effect.
#'   For a cubic cyclic effect, constructs the original \pkg{mgcv} cyclic
#'   spline basis and penalty matrix using the observed temporal range.
#'   \item Standardises continuous covariates using their observed means and
#'   standard deviations.
#'   \item Calls [create_insect_nimble_model()],
#'   [create_insect_nimble_data()], [create_insect_nimble_constants()], and
#'   [create_insect_nimble_initial_values()] (once per chain, each with its
#'   own seed).
#'   \item Sources the generated model code, builds the NIMBLE model (using
#'   the first chain's initial values), and compiles it.
#'   \item If `run_mcmc = TRUE`, configures and compiles the MCMC sampler and
#'   runs the requested iterations with chain-specific initial values and
#'   seeds.
#' }
#'
#' @return A named list containing the components used to construct and fit
#'   the model, including:
#' \describe{
#'   \item{`model_file`}{Information returned by
#'   [create_insect_nimble_model()].}
#'   \item{`model_code`}{NIMBLE model code sourced from the generated file.}
#'   \item{`data`}{Prepared response and predictor data supplied to NIMBLE.}
#'   \item{`constants`}{Model constants.}
#'   \item{`initial_values`}{List of length `nchains`; element `i` holds the
#'   initial values used for chain `i`.}
#'   \item{`model`, `compiled_model`}{NIMBLE model objects (returned when
#'   `run_mcmc = FALSE`).}
#'   \item{`samples`}{When `run_mcmc = TRUE`, the output of
#'   [nimble::runMCMC()].}
#'   \item{`seed`}{The `seed` argument as supplied.}
#'   \item{`chain_seeds`}{Integer vector of the seed used for each chain.}
#'   \item{`package_versions`}{Named character vector of versions of the
#'   packages used (`NA` if a package is not installed).}
#'   \item{`session`}{R version, platform, RNG kind and run timestamp.}
#'   \item{Additional metadata}{Response, distribution, temporal, covariate
#'   and spatial settings.}
#' }
#'
#' @seealso
#' [create_insect_nimble_model()],
#' [create_insect_nimble_data()],
#' [create_insect_nimble_constants()],
#' [create_insect_nimble_initial_values()]
#'
#' @importFrom stats na.omit
#' @export
run_insect_nimble_model <- function(
  modelling_data,
  modelCovs,
  covariate_types = NULL,
  share_fixed_effects = FALSE,
  share_hyperparameters = FALSE,

  # Responses
  response = "both",
  biomass_col = "sum_wet_weight",
  diversity_col = "shannon_div",

  # Response distributions
  biomass_distribution = "normal",
  diversity_distribution = "normal",

  # Replicate structure
  replicate_col = "sampling_name",

  # Temporal structure
  temporal_col = c("start_date_julian", "end_date_julian"),

  # Random effects
  random_year = TRUE,
  random_locality = TRUE,
  random_replicate = TRUE,

  # Temporal effect
  temporal_effect = c("cubic_cyclic", "none", "linear"),
  temporal_shared = TRUE,
  temporal_basis_k = 10,

  # Spatial effect: "none", "coordinates", "autocorrelation"
  # If spatial_shared = TRUE:
  #   - coordinates: one shared spatial coordinate effect
  #   - autocorrelation: one shared latent spatial process
  # If FALSE: biomass and diversity receive separate spatial effects
  spatial_effect = c("none", "coordinates", "autocorrelation"),
  spatial_shared = TRUE,
  longitude_col = "longitude",
  latitude_col = "latitude",

  # Spatial priors. Range is in metres because coordinates are
  # transformed to EPSG:25833 for the distance calculation.
  spatial_range_prior = c(1000, 500000),
  spatial_sigma_prior = 1,

  # Output
  output_dir = tempdir(),

  # MCMC
  seed = 123,
  run_mcmc = TRUE,
  run_waic = TRUE,
  niter = 20000,
  nburnin = 5000,
  thin = 10,
  nchains = 2
) {

  # ==========================================================
  # Load required package
  # ==========================================================

  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required.", call. = FALSE)
  }

  # ==========================================================
  # Validate chains and seeds
  # ==========================================================

  if (!is.numeric(nchains) ||
      length(nchains) != 1L ||
      !is.finite(nchains) ||
      nchains < 1 ||
      nchains != as.integer(nchains)) {
    stop("'nchains' must be a single positive integer.", call. = FALSE)
  }

  nchains <- as.integer(nchains)

  if (!is.numeric(seed) ||
      !(length(seed) %in% c(1L, nchains)) ||
      any(!is.finite(seed)) ||
      any(seed != round(seed))) {
    stop(
      "'seed' must be a single integer or an integer vector ",
      "of length 'nchains' (", nchains, ").",
      call. = FALSE
    )
  }

  if (length(seed) == 1L) {
    chain_seeds <- as.integer(seed) + seq_len(nchains) - 1L
  } else {
    chain_seeds <- as.integer(seed)
  }

  if (anyDuplicated(chain_seeds)) {
    stop("Each chain must have a distinct seed.", call. = FALSE)
  }

  names(chain_seeds) <- paste0("chain", seq_len(nchains))

  # ==========================================================
  # Validate response
  # ==========================================================

  response <- match.arg(
    response,
    choices = c("biomass", "diversity", "both")
  )

  # ==========================================================
  # Validate distributions
  # ==========================================================

  biomass_distribution <- match.arg(
    biomass_distribution,
    choices = c("normal", "lognormal")
  )

  diversity_distribution <- match.arg(
    diversity_distribution,
    choices = c("normal", "poisson", "negative_binomial")
  )

  # ==========================================================
  # Validate temporal sharing
  # ==========================================================

  if (!is.logical(temporal_shared) ||
      length(temporal_shared) != 1L ||
      is.na(temporal_shared)) {
    stop("'temporal_shared' must be TRUE or FALSE.", call. = FALSE)
  }

  # ==========================================================
  # Validate spatial effect
  # ==========================================================

  spatial_effect <- match.arg(
    spatial_effect,
    choices = c("none", "coordinates", "autocorrelation")
  )

  # ==========================================================
  # Validate spatial priors
  # ==========================================================

  if (length(spatial_range_prior) != 2 ||
      any(!is.finite(spatial_range_prior)) ||
      spatial_range_prior[1] <= 0 ||
      spatial_range_prior[2] <= spatial_range_prior[1]) {
    stop(
      "'spatial_range_prior' must contain two finite ",
      "positive values with upper > lower.",
      call. = FALSE
    )
  }

  if (length(spatial_sigma_prior) != 1 ||
      !is.finite(spatial_sigma_prior) ||
      spatial_sigma_prior <= 0) {
    stop(
      "'spatial_sigma_prior' must be a single positive finite value.",
      call. = FALSE
    )
  }

  # ==========================================================
  # Validate modelling data
  # ==========================================================

  if (!is.data.frame(modelling_data)) {
    stop("'modelling_data' must be a data.frame.", call. = FALSE)
  }

  if (nrow(modelling_data) == 0) {
    stop("'modelling_data' contains no rows.", call. = FALSE)
  }

  # ==========================================================
  # Validate model covariates
  # ==========================================================

  missing_covariates <- setdiff(modelCovs, names(modelling_data))

  if (length(missing_covariates) > 0) {
    stop(
      "The following model covariates were not found in ",
      "'modelling_data': ",
      paste(missing_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  # ==========================================================
  # Validate covariate types
  # ==========================================================

  if (is.null(covariate_types)) {
    covariate_types <- setNames(
      rep("continuous", length(modelCovs)),
      modelCovs
    )
  }

  missing_types <- setdiff(modelCovs, names(covariate_types))

  if (length(missing_types) > 0) {
    stop(
      "No covariate type supplied for: ",
      paste(missing_types, collapse = ", "),
      call. = FALSE
    )
  }

  invalid_types <- setdiff(
    unique(covariate_types[modelCovs]),
    c("continuous", "categorical")
  )

  if (length(invalid_types) > 0) {
    stop(
      "Invalid covariate type(s): ",
      paste(invalid_types, collapse = ", "),
      ". Allowed types are 'continuous' and 'categorical'.",
      call. = FALSE
    )
  }

  # ==========================================================
  # Validate spatial coordinates (only when requested)
  # ==========================================================

  if (spatial_effect != "none") {

    required_spatial_columns <- c(longitude_col, latitude_col)

    missing_spatial_columns <- setdiff(
      required_spatial_columns,
      names(modelling_data)
    )

    if (length(missing_spatial_columns) > 0) {
      stop(
        "The following spatial coordinate columns were ",
        "not found in 'modelling_data': ",
        paste(missing_spatial_columns, collapse = ", "),
        call. = FALSE
      )
    }

    if (!is.numeric(modelling_data[[longitude_col]]) ||
        !is.numeric(modelling_data[[latitude_col]])) {
      stop("Spatial coordinate columns must be numeric.", call. = FALSE)
    }

    if (any(!is.finite(modelling_data[[longitude_col]])) ||
        any(!is.finite(modelling_data[[latitude_col]]))) {
      stop(
        "Spatial coordinates contain NA, NaN, or Inf values.",
        call. = FALSE
      )
    }
  }

  # ==========================================================
  # Extract response variables
  # ==========================================================

  biomass <- NULL
  diversity <- NULL

  if (response %in% c("biomass", "both")) {

    if (!biomass_col %in% names(modelling_data)) {
      stop("Biomass column not found: ", biomass_col, call. = FALSE)
    }

    biomass <- modelling_data[[biomass_col]]
  }

  if (response %in% c("diversity", "both")) {

    if (is.null(diversity_col)) {
      stop(
        "diversity_col must be supplied when ",
        "response includes diversity.",
        call. = FALSE
      )
    }

    if (!diversity_col %in% names(modelling_data)) {
      stop("Diversity column not found: ", diversity_col, call. = FALSE)
    }

    diversity <- modelling_data[[diversity_col]]
  }

  # ==========================================================
  # Replicate number
  # ==========================================================

  replicates <- NULL

  if (!is.null(replicate_col)) {

    if (!replicate_col %in% names(modelling_data)) {
      stop("Replicate column not found: ", replicate_col, call. = FALSE)
    }

    required_replicate_cols <- c("locality", "year")

    missing_replicate_cols <- setdiff(
      required_replicate_cols,
      names(modelling_data)
    )

    if (length(missing_replicate_cols) > 0) {
      stop(
        "The following columns are required to construct ",
        "replicate IDs: ",
        paste(missing_replicate_cols, collapse = ", "),
        call. = FALSE
      )
    }

    replicates <- modelling_data |>
      dplyr::mutate(
        replicate_number = as.numeric(
          sub(".*-", "", .data[[replicate_col]])
        )
      ) |>
      dplyr::group_by(locality, year) |>
      dplyr::mutate(
        sampling_id = as.numeric(as.factor(replicate_number))
      ) |>
      dplyr::ungroup() |>
      dplyr::pull(sampling_id)
  }

  # ==========================================================
  # Temporal effect
  # ==========================================================

  temporal_data <- NULL

  # Retained for the NIMBLE constants if needed
  temporal_basis_penalized <- NULL
  temporal_penalty_eigenvalues <- NULL
  n_temporal_basis <- 0L

  # Original mgcv objects
  temporal_basis <- NULL
  temporal_penalty <- NULL

  temporal_effect <- match.arg(
    temporal_effect,
    choices = c("none", "linear", "cubic_cyclic")
  )

  if (!is.numeric(temporal_basis_k) ||
      length(temporal_basis_k) != 1L ||
      !is.finite(temporal_basis_k) ||
      temporal_basis_k < 4 ||
      temporal_basis_k != as.integer(temporal_basis_k)) {
    stop("'temporal_basis_k' must be a single integer >= 4.", call. = FALSE)
  }

  temporal_basis_k <- as.integer(temporal_basis_k)

  if (temporal_effect != "none") {

    if (is.null(temporal_col)) {
      stop(
        "temporal_col must be supplied when ",
        "temporal_effect is not 'none'.",
        call. = FALSE
      )
    }

    missing_temporal <- setdiff(temporal_col, names(modelling_data))

    if (length(missing_temporal) > 0) {
      stop(
        "Temporal columns not found: ",
        paste(missing_temporal, collapse = ", "),
        call. = FALSE
      )
    }

    temporal_data <- rowMeans(
      modelling_data[, temporal_col, drop = FALSE],
      na.rm = TRUE
    )

    if (any(!is.finite(temporal_data))) {
      stop(
        "Temporal data contains NA, NaN, or Inf values.",
        call. = FALSE
      )
    }

    # ------------------------------------------------------
    # Linear temporal effect
    # ------------------------------------------------------

    if (temporal_effect == "linear") {
      temporal_data <- as.numeric(scale(temporal_data))
    }

    # ------------------------------------------------------
    # Cubic cyclic temporal effect
    # ------------------------------------------------------

    if (temporal_effect == "cubic_cyclic") {

      # 1. Validate dependencies and inputs

      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop(
          "Package 'mgcv' is required for ",
          "temporal_effect = 'cubic_cyclic'.",
          call. = FALSE
        )
      }

      if (!is.numeric(temporal_data) ||
          length(temporal_data) < 3L ||
          any(!is.finite(temporal_data))) {
        stop(
          "'temporal_data' must be a numeric vector containing ",
          "at least three finite observations.",
          call. = FALSE
        )
      }

      # 2. Preserve the original temporal scale. Do not
      #    standardise: the original scale defines the cyclic
      #    boundary conditions.

      temporal_smooth_data <- data.frame(
        temporal = as.numeric(temporal_data)
      )

      temporal_boundary <- range(temporal_smooth_data$temporal)

      if (length(temporal_boundary) != 2L ||
          !all(is.finite(temporal_boundary)) ||
          diff(temporal_boundary) <= 0) {
        stop(
          "The temporal data must span a non-zero finite range ",
          "to construct a cyclic cubic spline.",
          call. = FALSE
        )
      }

      # 3. Construct the cyclic cubic spline using mgcv

      temporal_smooth_spec <- mgcv::s(
        temporal,
        bs = "cc",
        k = temporal_basis_k
      )

      temporal_knots <- list(temporal = temporal_boundary)

      temporal_smooth <- mgcv::smoothCon(
        temporal_smooth_spec,
        data = temporal_smooth_data,
        knots = temporal_knots,
        absorb.cons = FALSE,
        scale.penalty = FALSE
      )[[1L]]

      # 4. Extract the original basis and penalty. The
      #    eigen-decomposition and projection are performed by
      #    create_insect_nimble_model().

      temporal_basis <- temporal_smooth$X
      temporal_penalty <- temporal_smooth$S[[1L]]

      # 5. Validate the basis and penalty

      if (!is.matrix(temporal_basis) || !is.numeric(temporal_basis)) {
        stop(
          "Failed to construct a valid numeric mgcv temporal basis.",
          call. = FALSE
        )
      }

      if (!is.matrix(temporal_penalty) || !is.numeric(temporal_penalty)) {
        stop(
          "Failed to construct a valid numeric mgcv temporal penalty.",
          call. = FALSE
        )
      }

      if (nrow(temporal_basis) != length(temporal_data)) {
        stop(
          "The temporal basis row count does not match ",
          "the length of 'temporal_data'.",
          call. = FALSE
        )
      }

      if (nrow(temporal_penalty) != ncol(temporal_basis) ||
          ncol(temporal_penalty) != ncol(temporal_basis)) {
        stop(
          "The temporal penalty dimensions must match ",
          "the number of columns in the temporal basis.",
          call. = FALSE
        )
      }

      if (any(!is.finite(temporal_basis)) ||
          any(!is.finite(temporal_penalty))) {
        stop(
          "The temporal basis or penalty contains non-finite values.",
          call. = FALSE
        )
      }

      # Check penalty symmetry within numerical tolerance
      penalty_scale <- max(1, max(abs(temporal_penalty)))

      if (max(abs(temporal_penalty - t(temporal_penalty))) >
          sqrt(.Machine$double.eps) * penalty_scale) {
        stop(
          "The mgcv temporal penalty is not symmetric ",
          "within numerical tolerance.",
          call. = FALSE
        )
      }
    }
  }

  # ==========================================================
  # Identify and scale continuous covariates
  # ==========================================================

  continuous_covariates <- modelCovs[
    covariate_types[modelCovs] == "continuous"
  ]

  if (length(continuous_covariates) > 0) {

    scale_covariate <- function(x) {

      if (!is.numeric(x)) {
        stop("Continuous covariates must be numeric.", call. = FALSE)
      }

      mean_x <- mean(x, na.rm = TRUE)
      sd_x <- sd(x, na.rm = TRUE)

      if (!is.finite(sd_x) || sd_x == 0) {
        stop(
          "A continuous covariate has zero or undefined ",
          "standard deviation.",
          call. = FALSE
        )
      }

      as.numeric((x - mean_x) / sd_x)
    }

    modelling_data <- modelling_data |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(continuous_covariates),
          scale_covariate
        )
      )
  }

  # ==========================================================
  # Create NIMBLE model file
  # ==========================================================

  model_file <- create_insect_nimble_model(
    modelCovs = modelCovs,
    random_year = random_year,
    random_locality = random_locality,
    random_replicate = random_replicate,
    temporal_effect = temporal_effect,
    temporal_shared = temporal_shared,
    # Supply the original mgcv basis and penalty
    temporal_basis = temporal_basis,
    temporal_penalty = temporal_penalty,
    output_dir = output_dir,
    covariate_types = covariate_types,
    response = response,
    biomass_distribution = biomass_distribution,
    diversity_distribution = diversity_distribution,
    # Spatial structure
    spatial_effect = spatial_effect,
    spatial_shared = spatial_shared,
    spatial_range_prior = spatial_range_prior,
    spatial_sigma_prior = spatial_sigma_prior,
    share_fixed_effects = share_fixed_effects,
    share_hyperparameters = share_hyperparameters
  )

  if (temporal_effect == "cubic_cyclic") {

    temporal_basis_penalized <- model_file$temporal_basis_penalized
    temporal_penalty_eigenvalues <- model_file$temporal_penalty_eigenvalues
    n_temporal_basis <- model_file$n_temporal_basis

    if (is.null(temporal_basis_penalized)) {
      stop(
        "The model constructor did not return ",
        "'temporal_basis_penalized'.",
        call. = FALSE
      )
    }

    if (!is.matrix(temporal_basis_penalized) ||
        nrow(temporal_basis_penalized) != length(temporal_data) ||
        ncol(temporal_basis_penalized) != n_temporal_basis ||
        length(temporal_penalty_eigenvalues) != n_temporal_basis) {
      stop(
        "The projected temporal basis or its eigenvalues ",
        "have inconsistent dimensions.",
        call. = FALSE
      )
    }
  }

  # ==========================================================
  # Create NIMBLE data
  # ==========================================================

  data <- create_insect_nimble_data(
    modelling_data = modelling_data,
    modelCovs = modelCovs,
    covariate_types = covariate_types,
    response = response,
    biomass_col = biomass_col,
    diversity_col = diversity_col,
    biomass_distribution = biomass_distribution,
    diversity_distribution = diversity_distribution,
    temporal_data = temporal_data,
    temporal_effect = temporal_effect,
    # Spatial structure
    spatial_effect = spatial_effect,
    longitude_col = longitude_col,
    latitude_col = latitude_col
  )

  # ==========================================================
  # Create NIMBLE constants
  #
  # For autocorrelation, this calculates ONLY the fixed
  # distance matrix. The spatial covariance matrix is built
  # inside NIMBLE from estimated rho and sigma.
  # ==========================================================

  constants <- create_insect_nimble_constants(
    modelling_data = modelling_data,
    modelCovs = modelCovs,
    covariate_types = covariate_types,
    replicate_col = replicate_col,
    random_year = random_year,
    random_locality = random_locality,
    random_replicate = random_replicate,
    # Spatial structure
    spatial_effect = spatial_effect,
    longitude_col = longitude_col,
    latitude_col = latitude_col
  )

  if (temporal_effect == "cubic_cyclic") {
    constants <- c(
      constants,
      list(
        n_temporal_basis = n_temporal_basis,
        temporal_basis_penalized = temporal_basis_penalized,
        temporal_penalty_eigenvalues = temporal_penalty_eigenvalues
      )
    )
  }

  # ==========================================================
  # Create initial values: one set per chain, each generated
  # with its own seed
  # ==========================================================

  initial_values <- lapply(
    chain_seeds,
    function(chain_seed) {
      create_insect_nimble_initial_values(
        modelling_data = modelling_data,
        modelCovs = modelCovs,
        covariate_types = covariate_types,
        share_fixed_effects = share_fixed_effects,
        share_hyperparameters = share_hyperparameters,
        response = response,
        biomass_col = biomass_col,
        diversity_col = diversity_col,
        biomass_distribution = biomass_distribution,
        diversity_distribution = diversity_distribution,
        random_year = random_year,
        random_locality = random_locality,
        random_replicate = random_replicate,
        temporal_data = temporal_data,
        temporal_basis_penalized = temporal_basis_penalized,
        temporal_effect = temporal_effect,
        temporal_shared = if (temporal_effect == "none") FALSE else temporal_shared,
        seed = chain_seed,
        # Spatial structure
        spatial_effect = spatial_effect,
        spatial_shared = spatial_shared,
        longitude_col = longitude_col,
        latitude_col = latitude_col,
        spatial_range_prior = spatial_range_prior,
        spatial_sigma_prior = spatial_sigma_prior
      )
    }
  )

  names(initial_values) <- names(chain_seeds)

  # Warn if the generator ignored the seed and produced identical
  # starting points for different chains
  if (nchains > 1L &&
      all(vapply(
        initial_values[-1L],
        identical,
        logical(1),
        initial_values[[1L]]
      ))) {
    warning(
      "All chains received identical initial values. Check that ",
      "create_insect_nimble_initial_values() uses its 'seed' argument.",
      call. = FALSE
    )
  }

  # ==========================================================
  # Reproducibility information
  # ==========================================================

  version_packages <- c("nimble", "mgcv", "dplyr", "coda", "sf")

  this_package <- tryCatch(
    utils::packageName(environment(sys.function())),
    error = function(e) NULL
  )

  if (!is.null(this_package)) {
    version_packages <- unique(c(this_package, version_packages))
  }

  package_versions <- vapply(
    version_packages,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  )

  session <- list(
    r_version = R.version.string,
    platform = R.version$platform,
    rng_kind = RNGkind(),
    run_started = Sys.time()
  )

  # ==========================================================
  # Read generated NIMBLE model
  # ==========================================================

  if (!file.exists(model_file$file)) {
    stop(
      "NIMBLE model file was not found:\n",
      model_file$file,
      call. = FALSE
    )
  }

  source(model_file$file, local = TRUE)

  model_code <- biomass_diversity_model

  # ==========================================================
  # Build NIMBLE model (with the first chain's initial values;
  # each chain's own values are supplied to runMCMC below)
  # ==========================================================

  nimble_model <- nimble::nimbleModel(
    code = model_code,
    constants = constants,
    data = data,
    inits = initial_values[[1L]],
    check = TRUE,
    calculate = FALSE
  )

  nimble_model$initializeInfo()

  # ==========================================================
  # Compile NIMBLE model
  # ==========================================================

  compiled_model <- nimble::compileNimble(nimble_model)

  # ==========================================================
  # Stop here if MCMC is not requested
  # ==========================================================

  if (!run_mcmc) {
    return(
      list(
        model_file = model_file,
        model_code = model_code,
        data = data,
        constants = constants,
        initial_values = initial_values,
        model = nimble_model,
        compiled_model = compiled_model,
        response = response,
        biomass_distribution = biomass_distribution,
        diversity_distribution = diversity_distribution,
        temporal_effect = temporal_effect,
        temporal_basis_k = temporal_basis_k,
        temporal_basis_penalized = temporal_basis_penalized,
        temporal_penalty_eigenvalues = temporal_penalty_eigenvalues,
        spatial_effect = spatial_effect,
        spatial_shared = spatial_shared,
        # Reproducibility
        seed = seed,
        chain_seeds = chain_seeds,
        nchains = nchains,
        package_versions = package_versions,
        session = session
      )
    )
  }

  # ==========================================================
  # Configure, build and compile MCMC
  # ==========================================================

  mcmc_conf <- nimble::configureMCMC(
    nimble_model,
    enableWAIC = run_waic
  )

  mcmc <- nimble::buildMCMC(mcmc_conf)

  compiled_mcmc <- nimble::compileNimble(
    mcmc,
    project = nimble_model
  )

  # ==========================================================
  # Run MCMC
  #
  # inits:   list of length nchains -> chain i starts from
  #          initial_values[[i]]
  # setSeed: numeric vector of length nchains -> R's seed is
  #          set to chain_seeds[i] at the start of chain i
  # ==========================================================

  samples <- nimble::runMCMC(
    compiled_mcmc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    inits = unname(initial_values),
    setSeed = unname(chain_seeds),
    samplesAsCodaMCMC = TRUE,
    summary = TRUE,
    WAIC = run_waic
  )

  # ==========================================================
  # Return everything
  # ==========================================================
return(
  list(
    model_file = model_file,
    model_code = model_code,
    data = data,
    constants = constants,
    initial_values = initial_values,
    samples = samples,
    response = response,
    biomass_distribution = biomass_distribution,
    diversity_distribution = diversity_distribution,
    temporal_effect = temporal_effect,
    temporal_shared = temporal_shared,
    temporal_basis_k = temporal_basis_k,
    temporal_basis_penalized = temporal_basis_penalized,
    temporal_penalty_eigenvalues = temporal_penalty_eigenvalues,
    modelCovs = modelCovs,
    covariate_types = covariate_types,
    # Spatial model information
    spatial_effect = spatial_effect,
    spatial_shared = spatial_shared,
    longitude_col = longitude_col,
    latitude_col = latitude_col,
    spatial_range_prior = spatial_range_prior,
    spatial_sigma_prior = spatial_sigma_prior,
    # MCMC settings
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    # Reproducibility
    seed = seed,
    chain_seeds = chain_seeds,
    package_versions = package_versions,
    session = session
  )
)
}
