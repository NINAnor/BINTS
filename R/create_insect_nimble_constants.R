# ==============================================================================
# Create NIMBLE model constants
#
# Creates the constants required by the generated NIMBLE model.
#
# Continuous covariates:
#   - counted in p
#   - supplied through the X matrix in the data object
#
# Categorical covariates:
#   - converted to integer IDs
#   - supplied as <covariate>_id
#   - number of levels supplied as n_<covariate>
#
# Optional random effects:
#   - year
#   - locality
#   - replicate
#
# Optional spatial effects:
#   - none
#   - coordinates
#   - autocorrelation
#
# Spatial autocorrelation:
#   - spatial locations are defined by longitude/latitude
#   - coordinates are transformed to EPSG:25833
#   - distances are calculated in metres
#   - only the distance matrix is supplied to NIMBLE
#   - the spatial covariance matrix is constructed inside NIMBLE
#   - spatial range and spatial variance are estimated by NIMBLE
#
# IMPORTANT:
#
#   spatial_range and spatial_sigma are deliberately NOT calculated here.
#
#   They are model parameters and are estimated by NIMBLE.
#
# TEMPORAL EFFECTS:
#
#   Temporal basis functions and their penalty matrix are deliberately NOT
#   constructed here.
#
#   For temporal_effect = "cubic_cyclic", the mgcv basis and penalty are
#   transformed by create_insect_nimble_model() into the positive-eigenvalue
#   basis:
#
#       B_penalized = B %*% V_positive
#
#   with priors:
#
#       beta_k ~ Normal(
#           0,
#           [tau_temporal * lambda_k]^-1
#       )
#
#   where lambda_k are the positive eigenvalues of the penalty matrix.
#
# ==============================================================================


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