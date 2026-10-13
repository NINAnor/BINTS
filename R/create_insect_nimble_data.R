# ==============================================================================
# Create NIMBLE model data
#
# Prepares the data object required by the generated NIMBLE model.
#
# Response options:
#   - "biomass"
#   - "diversity"
#   - "both"
#
# Biomass distributions:
#   - "normal"
#   - "lognormal"
#
# Diversity distributions:
#   - "poisson"
#   - "negative_binomial"
#
# Continuous covariates:
#   -> supplied through X
#
# Categorical covariates:
#   -> represented through integer IDs in the constants object
#
# Temporal effects:
#   - "none"
#   - "linear"
#   - "cubic_cyclic"
#
# Spatial effects:
#   - "none"
#   - "coordinates"
#   - "autocorrelation"
#
# IMPORTANT:
#   The cubic cyclic spline basis and penalty matrix are NOT stored here.
#   They are supplied to the model-building function and transformed there
#   into the positive-eigenvalue basis used by the NIMBLE model.
#
# ==============================================================================


create_insect_nimble_data <- function(
    
    modelling_data,
    modelCovs,
    covariate_types = NULL,
    
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
    # Temporal effect
    # --------------------------------------------------------------------------
    
    temporal_data = NULL,
    temporal_effect = c(
        "none",
        "linear",
        "cubic_cyclic"
    ),
    
    # --------------------------------------------------------------------------
    # Spatial effect
    # --------------------------------------------------------------------------
    
    spatial_effect = c(
        "none",
        "coordinates",
        "autocorrelation"
    ),
    
    longitude_col = "longitude",
    latitude_col = "latitude"
    
) {
    
    # ==========================================================================
    # 1. Validate response
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
    # 2. Validate response distributions
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
    # 3. Validate temporal effect
    # ==========================================================================
    
    temporal_effect <- match.arg(
        temporal_effect,
        choices = c(
            "none",
            "linear",
            "cubic_cyclic"
        )
    )
    
    
    # ==========================================================================
    # 4. Validate spatial effect
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
    # 5. Determine included responses
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
    # 6. Validate modelling data
    # ==========================================================================
    
    if (!is.data.frame(modelling_data)) {
        
        stop(
            "modelling_data must be a data frame.",
            call. = FALSE
        )
        
    }
    
    N <- nrow(modelling_data)
    
    if (N == 0) {
        
        stop(
            "modelling_data contains no observations.",
            call. = FALSE
        )
        
    }
    
    
    # ==========================================================================
    # 7. Validate model covariates
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
    # 8. Validate covariate types
    # ==========================================================================
    
    if (length(modelCovs) > 0) {
        
        # ----------------------------------------------------------------------
        # Default: treat all covariates as continuous
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
            # Check that every model covariate has a type
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
    # 9. Identify continuous covariates
    # ==========================================================================
    
    continuous_covariates <- modelCovs[
        covariate_types[modelCovs] == "continuous"
    ]
    
    
    # ==========================================================================
    # 10. Identify categorical covariates
    # ==========================================================================
    
    categorical_covariates <- modelCovs[
        covariate_types[modelCovs] == "categorical"
    ]
    
    
    # ==========================================================================
    # 11. Initialise NIMBLE data object
    # ==========================================================================
    
    nimble_data <- list()
    
    
    # ==========================================================================
    # 12. Biomass response
    # ==========================================================================
    
    if (include_biomass) {
        
        # ----------------------------------------------------------------------
        # Check that biomass column was supplied
        # ----------------------------------------------------------------------
        
        if (is.null(biomass_col)) {
            
            stop(
                "biomass_col must be supplied when ",
                "response includes biomass.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check that biomass column exists
        # ----------------------------------------------------------------------
        
        if (!biomass_col %in% names(modelling_data)) {
            
            stop(
                "Biomass column not found in modelling_data: ",
                biomass_col,
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Extract biomass
        # ----------------------------------------------------------------------
        
        biomass <- modelling_data[[biomass_col]]
        
        # ----------------------------------------------------------------------
        # Numeric validation
        # ----------------------------------------------------------------------
        
        if (!is.numeric(biomass)) {
            
            stop(
                "Biomass response must be numeric.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Missing-value validation
        # ----------------------------------------------------------------------
        
        if (anyNA(biomass)) {
            
            stop(
                "Missing values found in biomass response: ",
                biomass_col,
                ". Remove or handle missing observations before fitting the model.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Finite-value validation
        # ----------------------------------------------------------------------
        
        if (any(!is.finite(biomass))) {
            
            stop(
                "Non-finite values found in biomass response: ",
                biomass_col,
                ".",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Lognormal validation
        # ----------------------------------------------------------------------
        
        if (biomass_distribution == "lognormal") {
            
            if (any(biomass <= 0)) {
                
                stop(
                    "Log-normal biomass requires strictly positive values. ",
                    "The biomass response contains ",
                    sum(biomass <= 0),
                    " value(s) <= 0.",
                    call. = FALSE
                )
                
            }
            
        }
        
        # ----------------------------------------------------------------------
        # Add biomass to NIMBLE data
        # ----------------------------------------------------------------------
        
        nimble_data[["biomass"]] <- as.numeric(biomass)
        
    }
    
    
    # ==========================================================================
    # 13. Diversity response
    # ==========================================================================
    
    if (include_diversity) {
        
        # ----------------------------------------------------------------------
        # Check that diversity column was supplied
        # ----------------------------------------------------------------------
        
        if (is.null(diversity_col)) {
            
            stop(
                "diversity_col must be supplied when ",
                "response includes diversity.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check that diversity column exists
        # ----------------------------------------------------------------------
        
        if (!diversity_col %in% names(modelling_data)) {
            
            stop(
                "Diversity column not found in modelling_data: ",
                diversity_col,
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Extract diversity
        # ----------------------------------------------------------------------
        
        diversity <- modelling_data[[diversity_col]]
        
        # ----------------------------------------------------------------------
        # Numeric validation
        # ----------------------------------------------------------------------
        
        if (!is.numeric(diversity)) {
            
            stop(
                "Diversity response must be numeric.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Missing-value validation
        # ----------------------------------------------------------------------
        
        if (anyNA(diversity)) {
            
            stop(
                "Missing values found in diversity response: ",
                diversity_col,
                ". Remove or handle missing observations before fitting the model.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Finite-value validation
        # ----------------------------------------------------------------------
        
        if (any(!is.finite(diversity))) {
            
            stop(
                "Non-finite values found in diversity response: ",
                diversity_col,
                ".",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Count-distribution validation
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
                    "The ",
                    diversity_distribution,
                    " distribution requires non-negative values. ",
                    "The diversity response contains ",
                    sum(diversity < 0),
                    " negative value(s).",
                    call. = FALSE
                )
                
            }
            
            if (any(diversity != floor(diversity))) {
                
                stop(
                    "The ",
                    diversity_distribution,
                    " distribution requires integer-valued observations. ",
                    "The diversity response contains non-integer values.",
                    call. = FALSE
                )
                
            }
            
        }
        
        # ----------------------------------------------------------------------
        # Add diversity to NIMBLE data
        # ----------------------------------------------------------------------
        
        nimble_data[["diversity"]] <- as.numeric(diversity)
        
    }
    
    
    # ==========================================================================
    # 14. Continuous covariates
    #
    # Only continuous covariates enter the X matrix.
    #
    # Categorical covariates are handled through integer IDs generated by
    # create_insect_nimble_constants().
    # ==========================================================================
    
    if (length(continuous_covariates) > 0) {
        
        # ----------------------------------------------------------------------
        # Check that all continuous covariates exist
        # ----------------------------------------------------------------------
        
        missing_covariates <- setdiff(
            continuous_covariates,
            names(modelling_data)
        )
        
        if (length(missing_covariates) > 0) {
            
            stop(
                "Continuous covariates not found in modelling_data: ",
                paste(
                    missing_covariates,
                    collapse = ", "
                ),
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Validate continuous covariates
        # ----------------------------------------------------------------------
        
        for (cov in continuous_covariates) {
            
            if (!is.numeric(modelling_data[[cov]])) {
                
                stop(
                    "Continuous covariate '",
                    cov,
                    "' must be numeric.",
                    call. = FALSE
                )
                
            }
            
            if (anyNA(modelling_data[[cov]])) {
                
                stop(
                    "Missing values found in continuous covariate: ",
                    cov,
                    ".",
                    call. = FALSE
                )
                
            }
            
            if (any(!is.finite(modelling_data[[cov]]))) {
                
                stop(
                    "Non-finite values found in continuous covariate: ",
                    cov,
                    ".",
                    call. = FALSE
                )
                
            }
            
        }
        
        # ----------------------------------------------------------------------
        # Construct design matrix
        # ----------------------------------------------------------------------
        
        X <- modelling_data[
            continuous_covariates
        ]
        
        X <- as.matrix(X)
        
        storage.mode(X) <- "double"
        
        # ----------------------------------------------------------------------
        # Verify dimensions
        # ----------------------------------------------------------------------
        
        if (
            nrow(X) != N ||
            ncol(X) != length(continuous_covariates)
        ) {
            
            stop(
                "Dimensions of X do not match the expected number ",
                "of observations and continuous covariates.",
                call. = FALSE
            )
            
        }
        
        nimble_data[["X"]] <- X
        
    } else {
        
        # ----------------------------------------------------------------------
        # No continuous covariates
        # ----------------------------------------------------------------------
        
        nimble_data[["X"]] <- matrix(
            numeric(0),
            nrow = N,
            ncol = 0
        )
        
    }
    
    
    # ==========================================================================
    # 15. Temporal variable
    #
    # For all temporal models:
    #
    #   temporal_effect = "linear"
    #       temporal is used directly as a linear predictor.
    #
    #   temporal_effect = "cubic_cyclic"
    #       temporal is used to construct/evaluate the cyclic spline basis.
    #       The basis and penalty are handled by the model-building function.
    #
    #   temporal_effect = "none"
    #       no temporal variable is added.
    #
    # ==========================================================================
    
    if (temporal_effect != "none") {
        
        # ----------------------------------------------------------------------
        # Check that temporal data were supplied
        # ----------------------------------------------------------------------
        
        if (is.null(temporal_data)) {
            
            stop(
                "temporal_data must be supplied when ",
                "temporal_effect = '",
                temporal_effect,
                "'.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check length
        # ----------------------------------------------------------------------
        
        if (length(temporal_data) != N) {
            
            stop(
                "temporal_data must have the same number of observations ",
                "as modelling_data. ",
                "Expected: ",
                N,
                "; supplied: ",
                length(temporal_data),
                ".",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check numeric
        # ----------------------------------------------------------------------
        
        if (!is.numeric(temporal_data)) {
            
            stop(
                "temporal_data must be numeric.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check missing values
        # ----------------------------------------------------------------------
        
        if (anyNA(temporal_data)) {
            
            stop(
                "Missing values found in temporal_data.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check finite values
        # ----------------------------------------------------------------------
        
        if (any(!is.finite(temporal_data))) {
            
            stop(
                "Non-finite values found in temporal_data.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Add temporal variable
        # ----------------------------------------------------------------------
        
        nimble_data[["temporal"]] <- as.numeric(
            temporal_data
        )
        
    }
    
    
    # ==========================================================================
    # 16. Spatial coordinates
    #
    # Coordinates are required when a spatial effect is requested.
    #
    # "coordinates":
    #   Longitude and latitude are included as standardised fixed effects.
    #
    # "autocorrelation":
    #   Longitude and latitude are validated here, but the original geographic
    #   coordinates remain in modelling_data.
    #
    #   create_insect_nimble_constants() is responsible for:
    #       1. transforming coordinates to EPSG:25833;
    #       2. identifying unique spatial locations;
    #       3. constructing spatial_id;
    #       4. calculating distances in metres.
    #
    # ==========================================================================
    
    spatial_coordinate_center <- NULL
    spatial_coordinate_scale <- NULL
    
    if (
        spatial_effect %in%
        c(
            "coordinates",
            "autocorrelation"
        )
    ) {
        
        # ----------------------------------------------------------------------
        # Check longitude column
        # ----------------------------------------------------------------------
        
        if (!longitude_col %in% names(modelling_data)) {
            
            stop(
                "Longitude column not found in modelling_data: ",
                longitude_col,
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Check latitude column
        # ----------------------------------------------------------------------
        
        if (!latitude_col %in% names(modelling_data)) {
            
            stop(
                "Latitude column not found in modelling_data: ",
                latitude_col,
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Extract coordinates
        # ----------------------------------------------------------------------
        
        longitude <- modelling_data[[longitude_col]]
        latitude <- modelling_data[[latitude_col]]
        
        # ----------------------------------------------------------------------
        # Numeric validation
        # ----------------------------------------------------------------------
        
        if (!is.numeric(longitude)) {
            
            stop(
                "Longitude must be numeric.",
                call. = FALSE
            )
            
        }
        
        if (!is.numeric(latitude)) {
            
            stop(
                "Latitude must be numeric.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Missing-value validation
        # ----------------------------------------------------------------------
        
        if (anyNA(longitude)) {
            
            stop(
                "Missing values found in longitude.",
                call. = FALSE
            )
            
        }
        
        if (anyNA(latitude)) {
            
            stop(
                "Missing values found in latitude.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Finite-value validation
        # ----------------------------------------------------------------------
        
        if (any(!is.finite(longitude))) {
            
            stop(
                "Non-finite values found in longitude.",
                call. = FALSE
            )
            
        }
        
        if (any(!is.finite(latitude))) {
            
            stop(
                "Non-finite values found in latitude.",
                call. = FALSE
            )
            
        }
        
        # ----------------------------------------------------------------------
        # Geographic range validation
        # ----------------------------------------------------------------------
        
        if (any(longitude < -180 | longitude > 180)) {
            
            stop(
                "Longitude values must be between -180 and 180 degrees.",
                call. = FALSE
            )
            
        }
        
        if (any(latitude < -90 | latitude > 90)) {
            
            stop(
                "Latitude values must be between -90 and 90 degrees.",
                call. = FALSE
            )
            
        }
        
        
        # ----------------------------------------------------------------------
        # Coordinate-based spatial effect
        # ----------------------------------------------------------------------
        
        if (spatial_effect == "coordinates") {
            
            longitude_mean <- mean(
                longitude
            )
            
            latitude_mean <- mean(
                latitude
            )
            
            longitude_sd <- stats::sd(
                longitude
            )
            
            latitude_sd <- stats::sd(
                latitude
            )
            
            # ------------------------------------------------------------------
            # Check that coordinates have sufficient variation
            # ------------------------------------------------------------------
            
            if (
                !is.finite(longitude_sd) ||
                longitude_sd == 0
            ) {
                
                stop(
                    "Longitude must have non-zero variation for ",
                    "spatial_effect = 'coordinates'.",
                    call. = FALSE
                )
                
            }
            
            if (
                !is.finite(latitude_sd) ||
                latitude_sd == 0
            ) {
                
                stop(
                    "Latitude must have non-zero variation for ",
                    "spatial_effect = 'coordinates'.",
                    call. = FALSE
                )
                
            }
            
            # ------------------------------------------------------------------
            # Standardise coordinates
            # ------------------------------------------------------------------
            
            nimble_data[["longitude"]] <- as.numeric(
                (
                    longitude - longitude_mean
                ) / longitude_sd
            )
            
            nimble_data[["latitude"]] <- as.numeric(
                (
                    latitude - latitude_mean
                ) / latitude_sd
            )
            
            # ------------------------------------------------------------------
            # Store transformation information
            # ------------------------------------------------------------------
            
            spatial_coordinate_center <- c(
                longitude = longitude_mean,
                latitude = latitude_mean
            )
            
            spatial_coordinate_scale <- c(
                longitude = longitude_sd,
                latitude = latitude_sd
            )
            
        } else {
            
            # ------------------------------------------------------------------
            # Autocorrelation model
            #
            # Do NOT add standardised coordinates to NIMBLE data.
            #
            # The constants function uses the original coordinates in
            # modelling_data.
            # ------------------------------------------------------------------
            
            spatial_coordinate_center <- NULL
            spatial_coordinate_scale <- NULL
            
        }
        
    }
    
    
    # ==========================================================================
    # 17. Store model-data metadata
    # ==========================================================================
    
    nimble_data_attributes <- list(
        
        response = response,
        
        biomass_col = biomass_col,
        
        diversity_col = diversity_col,
        
        biomass_distribution = biomass_distribution,
        
        diversity_distribution = diversity_distribution,
        
        continuous_covariates = continuous_covariates,
        
        categorical_covariates = categorical_covariates,
        
        temporal_effect = temporal_effect,
        
        spatial_effect = spatial_effect,
        
        longitude_col = longitude_col,
        
        latitude_col = latitude_col,
        
        spatial_coordinate_center = spatial_coordinate_center,
        
        spatial_coordinate_scale = spatial_coordinate_scale
        
    )
    
    
    attr(
        nimble_data,
        "model_metadata"
    ) <- nimble_data_attributes
    
    
    # ==========================================================================
    # 18. Return NIMBLE data
    # ==========================================================================
    
    nimble_data
    
}