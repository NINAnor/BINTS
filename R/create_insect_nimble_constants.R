#' Create Constants for an Insect NIMBLE Model
#'
#' Prepare the constants required by a generated NIMBLE model for insect
#' biomass, diversity, or both responses. The function derives the number of
#' observations, the number of continuous covariates, categorical covariate
#' indices, optional grouping-level indices, and spatial constants from the
#' modelling data.
#'
#' Continuous covariates are represented in the model design matrix and
#' contribute to the fixed-effect dimension `p`. Categorical covariates are
#' converted to integer indices, with the number of levels supplied through
#' corresponding constants. Optional random effects can be indexed by year,
#' locality, and replicate.
#'
#' When spatial autocorrelation is requested, geographic coordinates are
#' transformed from longitude and latitude in WGS84 (EPSG:4326) to
#' ETRS89 / UTM zone 33N (EPSG:25833). Pairwise Euclidean distances between
#' unique projected locations are calculated in metres. The function supplies
#' the distance matrix and spatial indices to NIMBLE; the spatial covariance
#' matrix and its range and standard-deviation parameters are constructed or
#' estimated within the generated model.
#'
#' Temporal spline basis functions and penalty matrices are not constructed by
#' this function. They are prepared separately and passed to
#' [create_insect_nimble_model()] when a penalized cyclic temporal effect is
#' required.
#'
#' @param modelling_data A data frame containing the observations used to
#'   construct the model constants. Required columns depend on the selected
#'   covariates and random or spatial effects.
#'
#' @param modelCovs A character vector specifying the covariates included in
#'   the model. If omitted or `NULL`, no covariates are included. Covariate
#'   names must correspond to columns in `modelling_data`.
#'
#' @param covariate_types An optional named character vector identifying each
#'   covariate as `"continuous"` or `"categorical"`. Names must correspond to
#'   elements of `modelCovs`. If `NULL`, all specified covariates are treated
#'   as continuous.
#'
#' @param replicate_col A character string specifying the column in
#'   `modelling_data` containing replicate identifiers. Required when
#'   `random_replicate = TRUE`. The replicate number is extracted from the
#'   trailing numeric component following the final hyphen. Replicate indices
#'   are assigned within locality and year.
#'
#' @param random_year Logical. If `TRUE`, construct integer year indices and
#'   the number of year levels for a year-level random effect. Defaults to
#'   `FALSE`.
#'
#' @param random_locality Logical. If `TRUE`, construct integer locality
#'   indices and the number of locality levels for a locality-level random
#'   effect. Defaults to `FALSE`.
#'
#' @param random_replicate Logical. If `TRUE`, construct replicate indices
#'   within locality and year, together with the number of replicate levels.
#'   Requires `replicate_col` and the `locality` and `year` columns.
#'   Defaults to `FALSE`.
#'
#' @param spatial_effect A character string specifying the spatial structure.
#'   Supported values are `"none"`, `"coordinates"`, and `"autocorrelation"`.
#'   The default is `"none"`. Spatial constants are constructed only for
#'   `"autocorrelation"`; coordinate-based effects require coordinates to be
#'   supplied as model data.
#'
#' @param spatial_shared Logical. Specifies whether the latent spatial
#'   autocorrelation process is shared between biomass and diversity.
#'   If `TRUE`, a common zero vector is generated for the shared process.
#'   If `FALSE`, separate zero vectors are generated for biomass and diversity.
#'   This argument affects the constants only when
#'   `spatial_effect = "autocorrelation"`. Defaults to `TRUE`.
#'
#' @param longitude_col A character string naming the column containing
#'   longitude coordinates in decimal degrees. Defaults to `"longitude"`.
#'   Used when `spatial_effect = "autocorrelation"`.
#'
#' @param latitude_col A character string naming the column containing
#'   latitude coordinates in decimal degrees. Defaults to `"latitude"`.
#'   Used when `spatial_effect = "autocorrelation"`.
#'
#' @return A named list of constants suitable for use in the generated NIMBLE
#'   model. The list contains:
#' \describe{
#'   \item{N}{Number of observations in `modelling_data`.}
#'   \item{p}{Number of continuous covariates.}
#'   \item{<covariate>_id}{Integer indices for each categorical covariate,
#'     using a syntactically valid covariate name.}
#'   \item{n_<covariate>}{Number of levels for each categorical covariate.}
#'   \item{year_id}{Integer year indices, if `random_year = TRUE`.}
#'   \item{n_year}{Number of year levels, if `random_year = TRUE`.}
#'   \item{locality_id}{Integer locality indices, if
#'     `random_locality = TRUE`.}
#'   \item{n_locality}{Number of locality levels, if
#'     `random_locality = TRUE`.}
#'   \item{replicate_id}{Integer replicate indices within locality and year,
#'     if `random_replicate = TRUE`.}
#'   \item{n_replicate}{Number of replicate indices, if
#'     `random_replicate = TRUE`.}
#'   \item{n_spatial}{Number of unique spatial locations, if
#'     `spatial_effect = "autocorrelation"`.}
#'   \item{spatial_id}{Integer indices mapping observations to unique spatial
#'     locations, if spatial autocorrelation is requested.}
#'   \item{distance_matrix}{Symmetric pairwise distance matrix in metres
#'     between unique projected spatial locations.}
#'   \item{spatial_zero}{Zero vector for the shared latent spatial process,
#'     if `spatial_shared = TRUE`.}
#'   \item{spatial_zero_biomass, spatial_zero_diversity}{Separate zero vectors
#'     for the response-specific spatial processes, if
#'     `spatial_shared = FALSE`.}
#' }
#'
#' @details
#' The function validates the modelling data, covariate specifications,
#' grouping variables, and spatial coordinates. Categorical variables and
#' grouping variables are converted to factors before their integer indices
#' are created. Missing values in categorical variables, enabled grouping
#' variables, continuous covariates, or spatial coordinates result in errors.
#'
#' Replicate identifiers are extracted from the final hyphen-delimited
#' component of the selected replicate column. These replicate numbers are
#' converted into integer sampling indices separately within each
#' locality-year combination.
#'
#' For spatial autocorrelation, observations at the same projected location
#' share a spatial index. Coordinates are rounded to six decimal places in
#' projected units when identifying unique locations. The resulting distance
#' matrix is symmetrized and its diagonal is set to zero.
#'
#' The spatial range parameter, spatial standard deviation, and spatial
#' covariance matrix are not estimated or constructed by this function.
#' They are handled by the generated NIMBLE model.
#'
#' The returned list provides constants only. The response variables,
#' continuous covariate design matrix, coordinates for coordinate-based
#' effects, and any required temporal basis information must be supplied
#' separately in the model data or model specification.
#'
#' @examples
#' \dontrun{
#' # Prepare constants for a model with continuous and categorical covariates.
#' model_constants <- create_insect_nimble_constants(
#'   modelling_data = insect_data,
#'   modelCovs = c("temperature", "habitat"),
#'   covariate_types = c(
#'     temperature = "continuous",
#'     habitat = "categorical"
#'   )
#' )
#'
#' # Inspect the number of observations and covariates.
#' model_constants$N
#' model_constants$p
#'
#' # Prepare constants for a model with spatial autocorrelation.
#' spatial_constants <- create_insect_nimble_constants(
#'   modelling_data = insect_data,
#'   modelCovs = c("temperature", "habitat"),
#'   covariate_types = c(
#'     temperature = "continuous",
#'     habitat = "categorical"
#'   ),
#'   spatial_effect = "autocorrelation",
#'   longitude_col = "longitude",
#'   latitude_col = "latitude"
#' )
#'
#' spatial_constants$n_spatial
#' dim(spatial_constants$distance_matrix)
#' }
#'
#' @seealso [create_insect_nimble_model()]
#' @export


create_insect_nimble_constants <- function(
    
    modelling_data,
    modelCovs,
    covariate_types = NULL,
    
    replicate_col = NULL,
    
    random_year = FALSE,
    random_locality = FALSE,
    random_replicate = FALSE,
    
    spatial_effect = c(
        "none",
        "coordinates",
        "autocorrelation"
    ),
    
    spatial_shared = TRUE,
    
    longitude_col = "longitude",
    latitude_col = "latitude"
    
) {
    
    
    # ==========================================================================
    # 0. Validate spatial arguments
    # ==========================================================================
    
    spatial_effect <- match.arg(
        spatial_effect
    )
    
    if (
        !is.logical(spatial_shared) ||
        length(spatial_shared) != 1 ||
        is.na(spatial_shared)
    ) {
        
        stop(
            "'spatial_shared' must be TRUE or FALSE.",
            call. = FALSE
        )
        
    }
    
    
    # ==========================================================================
    # 1. Validate random-effect arguments
    # ==========================================================================
    
    if (
        !is.logical(random_year) ||
        length(random_year) != 1 ||
        is.na(random_year)
    ) {
        
        stop(
            "'random_year' must be TRUE or FALSE.",
            call. = FALSE
        )
        
    }
    
    if (
        !is.logical(random_locality) ||
        length(random_locality) != 1 ||
        is.na(random_locality)
    ) {
        
        stop(
            "'random_locality' must be TRUE or FALSE.",
            call. = FALSE
        )
        
    }
    
    if (
        !is.logical(random_replicate) ||
        length(random_replicate) != 1 ||
        is.na(random_replicate)
    ) {
        
        stop(
            "'random_replicate' must be TRUE or FALSE.",
            call. = FALSE
        )
        
    }
    
    
    # ==========================================================================
    # 2. Validate modelling data
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
    # 3. Validate model covariates
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
    # 4. Validate covariate types
    # ==========================================================================
    
    if (length(modelCovs) > 0) {
        
        # ----------------------------------------------------------------------
        # Default: all covariates are continuous
        # ----------------------------------------------------------------------
        
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
            
            # ------------------------------------------------------------------
            # Check valid types
            # ------------------------------------------------------------------
            
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
            
            # ------------------------------------------------------------------
            # Check every model covariate has a type
            # ------------------------------------------------------------------
            
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
    # 5. Identify continuous and categorical covariates
    # ==========================================================================
    
    continuous_covariates <- modelCovs[
        covariate_types[modelCovs] == "continuous"
    ]
    
    categorical_covariates <- modelCovs[
        covariate_types[modelCovs] == "categorical"
    ]
    
    
    # ==========================================================================
    # 6. Initialise constants
    # ==========================================================================
    
    constants <- list(
        N = N,
        p = length(
            continuous_covariates
        )
    )
    
    
    # ==========================================================================
    # 7. Categorical covariates
    #
    # Each categorical covariate receives:
    #
    #   <covariate>_id
    #   n_<covariate>
    #
    # Example:
    #
    #   region_name_id
    #   n_region_name
    #
    # These correspond to the categorical random effects generated by
    # create_insect_nimble_model().
    # ==========================================================================
    
    if (length(categorical_covariates) > 0) {
        
        for (cov in categorical_covariates) {
            
            # ------------------------------------------------------------------
            # Check that covariate exists
            # ------------------------------------------------------------------
            
            if (!cov %in% names(modelling_data)) {
                
                stop(
                    "Categorical covariate not found in modelling_data: ",
                    cov,
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Convert to factor
            # ------------------------------------------------------------------
            
            cov_factor <- factor(
                modelling_data[[cov]]
            )
            
            # ------------------------------------------------------------------
            # Check missing values
            # ------------------------------------------------------------------
            
            if (anyNA(cov_factor)) {
                
                stop(
                    "Missing values found in categorical covariate: ",
                    cov,
                    ". Categorical random-effect IDs cannot contain NA.",
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Check number of levels
            # ------------------------------------------------------------------
            
            if (nlevels(cov_factor) < 2) {
                
                warning(
                    "Categorical covariate '",
                    cov,
                    "' has fewer than two levels. ",
                    "The associated random effect cannot estimate ",
                    "between-level variation.",
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Convert name to a valid NIMBLE object name
            # ------------------------------------------------------------------
            
            cov_clean <- make.names(
                cov
            )
            
            # ------------------------------------------------------------------
            # Integer ID for each observation
            # ------------------------------------------------------------------
            
            constants[[
                paste0(
                    cov_clean,
                    "_id"
                )
            ]] <- as.integer(
                cov_factor
            )
            
            # ------------------------------------------------------------------
            # Number of levels
            # ------------------------------------------------------------------
            
            constants[[
                paste0(
                    "n_",
                    cov_clean
                )
            ]] <- nlevels(
                cov_factor
            )
            
        }
        
    }
    
    
    # ==========================================================================
    # 8. Year random effect
    # ==========================================================================
    
    if (random_year) {
        
        # ----------------------------------------------------------------------
        # Check that year exists
        # ----------------------------------------------------------------------
        
        if (!"year" %in% names(modelling_data)) {
            
            stop(
                "'year' column not found in modelling_data.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Convert year to factor
        # ----------------------------------------------------------------------
        
        year_factor <- factor(
            modelling_data$year
        )
        
        # ----------------------------------------------------------------------
        # Check missing values
        # ----------------------------------------------------------------------
        
        if (anyNA(year_factor)) {
            
            stop(
                "Missing values found in 'year'. ",
                "Year random-effect IDs cannot contain NA.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check number of levels
        # ----------------------------------------------------------------------
        
        if (nlevels(year_factor) < 2) {
            
            warning(
                "'year' contains only one level. ",
                "The year random effect will not estimate between-year variation.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Create year IDs
        # ----------------------------------------------------------------------
        
        constants[["year_id"]] <- as.integer(
            year_factor
        )
        
        constants[["n_year"]] <- nlevels(
            year_factor
        )
        
    }
    
    
    # ==========================================================================
    # 9. Locality random effect
    # ==========================================================================
    
    if (random_locality) {
        
        # ----------------------------------------------------------------------
        # Check that locality exists
        # ----------------------------------------------------------------------
        
        if (!"locality" %in% names(modelling_data)) {
            
            stop(
                "'locality' column not found in modelling_data.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Convert locality to factor
        # ----------------------------------------------------------------------
        
        locality_factor <- factor(
            modelling_data$locality
        )
        
        # ----------------------------------------------------------------------
        # Check missing values
        # ----------------------------------------------------------------------
        
        if (anyNA(locality_factor)) {
            
            stop(
                "Missing values found in 'locality'. ",
                "Locality random-effect IDs cannot contain NA.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check number of levels
        # ----------------------------------------------------------------------
        
        if (nlevels(locality_factor) < 2) {
            
            warning(
                "'locality' contains only one level. ",
                "The locality random effect will not estimate ",
                "between-locality variation.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Create locality IDs
        # ----------------------------------------------------------------------
        
        constants[["locality_id"]] <- as.integer(
            locality_factor
        )
        
        constants[["n_locality"]] <- nlevels(
            locality_factor
        )
        
    }
    
    
    # ==========================================================================
    # 10. Replicate random effect
    #
    # The replicate number is extracted from the end of replicate_col.
    #
    # Example:
    #
    #   Ost_2020_Skog_1_sampling-01 -> 1
    #   Ost_2020_Skog_1_sampling-02 -> 2
    #
    # Replicate IDs are then assigned within locality × year.
    #
    # This preserves the interpretation of replicate as the sampling occasion
    # within a locality and year.
    # ==========================================================================
    
    if (random_replicate) {
        
        # ----------------------------------------------------------------------
        # replicate_col must be supplied
        # ----------------------------------------------------------------------
        
        if (is.null(replicate_col)) {
            
            stop(
                "replicate_col must be supplied when ",
                "random_replicate = TRUE.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check replicate column exists
        # ----------------------------------------------------------------------
        
        if (!replicate_col %in% names(modelling_data)) {
            
            stop(
                "'",
                replicate_col,
                "' column not found in modelling_data.",
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
        # Check extraction
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
        
        replicate_data <- modelling_data |>
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
        # Check resulting IDs
        # ----------------------------------------------------------------------
        
        if (anyNA(replicate_data$sampling_id)) {
            
            stop(
                "Missing values were generated in replicate IDs.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Create replicate constants
        # ----------------------------------------------------------------------
        
        constants[["replicate_id"]] <-
            as.integer(
                replicate_data$sampling_id
            )
        
        constants[["n_replicate"]] <-
            max(
                replicate_data$sampling_id
            )
        
    }
    
    
    # ==========================================================================
    # 11. Spatial effect
    #
    # "none":
    #   No spatial constants are created.
    #
    # "coordinates":
    #   Longitude and latitude are supplied as data to the NIMBLE model.
    #   No additional spatial constants are required here.
    #
    # "autocorrelation":
    #
    #   A latent spatial process is indexed by unique geographic locations.
    #
    #   Raw longitude/latitude are transformed from WGS84 to EPSG:25833 so that
    #   distances are calculated in metres.
    #
    #   The distance matrix is constructed here, but the covariance matrix is NOT.
    #
    #   NIMBLE estimates:
    #
    #       rho_spatial
    #       sigma_spatial
    #
    #   and constructs:
    #
    #       Sigma_ij =
    #           sigma_spatial^2 *
    #           exp(-distance_ij / rho_spatial)
    #
    #   inside the NIMBLE model.
    #
    # ==========================================================================
    
    if (spatial_effect == "autocorrelation") {
        
        # ----------------------------------------------------------------------
        # Check coordinate columns
        # ----------------------------------------------------------------------
        
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
        
        # ----------------------------------------------------------------------
        # Extract coordinates
        # ----------------------------------------------------------------------
        
        longitude <- modelling_data[[longitude_col]]
        latitude <- modelling_data[[latitude_col]]
        
        # ----------------------------------------------------------------------
        # Validate coordinates
        # ----------------------------------------------------------------------
        
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
                "spatial_effect = 'autocorrelation'.",
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
        
        # ----------------------------------------------------------------------
        # Convert longitude/latitude to sf points
        # ----------------------------------------------------------------------
        
        spatial_points <- modelling_data |>
            dplyr::select(
                dplyr::all_of(
                    c(
                        longitude_col,
                        latitude_col
                    )
                )
            ) |>
            sf::st_as_sf(
                coords = c(
                    longitude_col,
                    latitude_col
                ),
                crs = 4326,
                remove = FALSE
            )
        
        # ----------------------------------------------------------------------
        # Transform to EPSG:25833
        # ----------------------------------------------------------------------
        
        spatial_points_projected <- spatial_points |>
            sf::st_transform(
                crs = 25833
            )
        
        # ----------------------------------------------------------------------
        # Extract projected coordinates
        # ----------------------------------------------------------------------
        
        xy <- sf::st_coordinates(
            spatial_points_projected
        )
        
        xy <- as.matrix(
            xy[, 1:2, drop = FALSE]
        )
        
        storage.mode(xy) <- "double"
        
        # ----------------------------------------------------------------------
        # Define unique spatial locations
        #
        # Multiple observations at the same geographic location share the same
        # latent spatial effect.
        # ----------------------------------------------------------------------
        
        spatial_key <- paste(
            round(
                xy[, 1],
                6
            ),
            round(
                xy[, 2],
                6
            ),
            sep = "_"
        )
        
        unique_spatial_keys <- unique(
            spatial_key
        )
        
        spatial_id <- match(
            spatial_key,
            unique_spatial_keys
        )
        
        # ----------------------------------------------------------------------
        # Extract unique spatial coordinates
        # ----------------------------------------------------------------------
        
        unique_xy <- xy[
            match(
                unique_spatial_keys,
                spatial_key
            ),
            ,
            drop = FALSE
        ]
        
        # ----------------------------------------------------------------------
        # Number of spatial locations
        # ----------------------------------------------------------------------
        
        n_spatial <- nrow(
            unique_xy
        )
        
        if (n_spatial < 1) {
            
            stop(
                "No unique spatial locations were identified.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Calculate pairwise distances in metres
        #
        # This is the ONLY spatial matrix calculated outside NIMBLE.
        # ----------------------------------------------------------------------
        
        if (n_spatial == 1) {
            
            distance_matrix <- matrix(
                0,
                nrow = 1,
                ncol = 1
            )
            
        } else {
            
            distance_matrix <- as.matrix(
                stats::dist(
                    unique_xy
                )
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Force exact symmetry and numerical stability
        # ----------------------------------------------------------------------
        
        distance_matrix <- (
            distance_matrix +
                t(distance_matrix)
        ) / 2
        
        diag(
            distance_matrix
        ) <- 0
        
        storage.mode(
            distance_matrix
        ) <- "double"
        
        # ----------------------------------------------------------------------
        # Validate distance matrix
        # ----------------------------------------------------------------------
        
        if (
            any(!is.finite(distance_matrix))
        ) {
            
            stop(
                "The spatial distance matrix contains non-finite values.",
                call. = FALSE
            )
            
        }
        
        if (
            any(distance_matrix < 0)
        ) {
            
            stop(
                "The spatial distance matrix contains negative distances.",
                call. = FALSE
            )
            
        }
        
        if (
            nrow(distance_matrix) != n_spatial ||
            ncol(distance_matrix) != n_spatial
        ) {
            
            stop(
                "The spatial distance matrix has incorrect dimensions.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Store spatial constants
        # ----------------------------------------------------------------------
        
        constants[["n_spatial"]] <-
            as.integer(
                n_spatial
            )
        
        constants[["spatial_id"]] <-
            as.integer(
                spatial_id
            )
        
        constants[["distance_matrix"]] <-
            distance_matrix
        
        # ----------------------------------------------------------------------
        # Zero vector for spatial MVN mean
        #
        # spatial_shared = TRUE:
        #   one common latent spatial process.
        #
        # spatial_shared = FALSE:
        #   separate latent processes for biomass and diversity.
        # ----------------------------------------------------------------------
        
        if (spatial_shared) {
            
            constants[["spatial_zero"]] <-
                rep(
                    0,
                    n_spatial
                )
            
        } else {
            
            constants[["spatial_zero_biomass"]] <-
                rep(
                    0,
                    n_spatial
                )
            
            constants[["spatial_zero_diversity"]] <-
                rep(
                    0,
                    n_spatial
                )
            
        }
        
    }
    
    
    # ==========================================================================
    # 12. Validate continuous covariates
    # ==========================================================================
    
    if (length(continuous_covariates) > 0) {
        
        for (cov in continuous_covariates) {
            
            # ------------------------------------------------------------------
            # Check that covariate exists
            # ------------------------------------------------------------------
            
            if (!cov %in% names(modelling_data)) {
                
                stop(
                    "Continuous covariate not found in modelling_data: ",
                    cov,
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Check numeric
            # ------------------------------------------------------------------
            
            if (!is.numeric(modelling_data[[cov]])) {
                
                stop(
                    "Continuous covariate '",
                    cov,
                    "' must be numeric.",
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Check missing values
            # ------------------------------------------------------------------
            
            if (anyNA(modelling_data[[cov]])) {
                
                stop(
                    "Missing values found in continuous covariate: ",
                    cov,
                    ". Missing predictor values are not allowed in the current model.",
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Check finite values
            # ------------------------------------------------------------------
            
            if (any(!is.finite(modelling_data[[cov]]))) {
                
                stop(
                    "Non-finite values found in continuous covariate: ",
                    cov,
                    ".",
                    call. = FALSE
                )
                
            }
            
        }
        
    }
    
    
    # ==========================================================================
    # 13. Validate dimensions of ID vectors
    # ==========================================================================
    
    
    # --------------------------------------------------------------------------
    # Categorical IDs
    # --------------------------------------------------------------------------
    
    if (length(categorical_covariates) > 0) {
        
        for (cov in categorical_covariates) {
            
            cov_clean <- make.names(
                cov
            )
            
            id_name <- paste0(
                cov_clean,
                "_id"
            )
            
            if (
                length(
                    constants[[id_name]]
                ) != N
            ) {
                
                stop(
                    "Length of '",
                    id_name,
                    "' does not equal N.",
                    call. = FALSE
                )
                
            }
            
            if (
                any(
                    constants[[id_name]] < 1 |
                        constants[[id_name]] >
                        constants[[
                            paste0(
                                "n_",
                                cov_clean
                            )
                        ]]
                )
            ) {
                
                stop(
                    "'",
                    id_name,
                    "' contains invalid category indices.",
                    call. = FALSE
                )
                
            }
            
        }
        
    }
    
    
    # --------------------------------------------------------------------------
    # Year ID
    # --------------------------------------------------------------------------
    
    if (random_year) {
        
        if (
            length(
                constants[["year_id"]]
            ) != N
        ) {
            
            stop(
                "Length of 'year_id' does not equal N.",
                call. = FALSE
            )
            
        }
        
        if (
            any(
                constants[["year_id"]] < 1 |
                    constants[["year_id"]] >
                    constants[["n_year"]]
            )
        ) {
            
            stop(
                "'year_id' contains invalid year indices.",
                call. = FALSE
            )
            
        }
        
    }
    
    
    # --------------------------------------------------------------------------
    # Locality ID
    # --------------------------------------------------------------------------
    
    if (random_locality) {
        
        if (
            length(
                constants[["locality_id"]]
            ) != N
        ) {
            
            stop(
                "Length of 'locality_id' does not equal N.",
                call. = FALSE
            )
            
        }
        
        if (
            any(
                constants[["locality_id"]] < 1 |
                    constants[["locality_id"]] >
                    constants[["n_locality"]]
            )
        ) {
            
            stop(
                "'locality_id' contains invalid locality indices.",
                call. = FALSE
            )
            
        }
        
    }
    
    
    # --------------------------------------------------------------------------
    # Replicate ID
    # --------------------------------------------------------------------------
    
    if (random_replicate) {
        
        if (
            length(
                constants[["replicate_id"]]
            ) != N
        ) {
            
            stop(
                "Length of 'replicate_id' does not equal N.",
                call. = FALSE
            )
            
        }
        
        if (
            any(
                constants[["replicate_id"]] < 1 |
                    constants[["replicate_id"]] >
                    constants[["n_replicate"]]
            )
        ) {
            
            stop(
                "'replicate_id' contains invalid replicate indices.",
                call. = FALSE
            )
            
        }
        
    }
    
    
    # --------------------------------------------------------------------------
    # Spatial ID
    # --------------------------------------------------------------------------
    
    if (spatial_effect == "autocorrelation") {
        
        if (
            length(
                constants[["spatial_id"]]
            ) != N
        ) {
            
            stop(
                "Length of 'spatial_id' does not equal N.",
                call. = FALSE
            )
            
        }
        
        if (
            any(
                constants[["spatial_id"]] < 1 |
                    constants[["spatial_id"]] >
                    constants[["n_spatial"]]
            )
        ) {
            
            stop(
                "'spatial_id' contains invalid spatial indices.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Validate distance matrix dimensions
        # ----------------------------------------------------------------------
        
        if (
            nrow(
                constants[["distance_matrix"]]
            ) != constants[["n_spatial"]] ||
            ncol(
                constants[["distance_matrix"]]
            ) != constants[["n_spatial"]]
        ) {
            
            stop(
                "'distance_matrix' dimensions do not match 'n_spatial'.",
                call. = FALSE
            )
            
        }
        
    }
    
    
    # ==========================================================================
    # 14. Return constants
    # ==========================================================================
    
    constants
    
}