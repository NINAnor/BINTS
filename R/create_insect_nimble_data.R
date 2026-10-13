#' Prepare Data for an Insect NIMBLE Model
#'
#' Construct the data list required by a generated NIMBLE model for insect
#' biomass, diversity, or both responses. The function extracts and validates
#' response variables, prepares the continuous covariate design matrix,
#' validates categorical covariates, and optionally prepares temporal data
#' and spatial coordinate predictors.
#'
#' Biomass can be modelled using a normal or log-normal observation
#' distribution, while diversity can be modelled using a Poisson or
#' negative-binomial observation distribution. The selected responses are
#' returned under the names `biomass` and `diversity`, respectively.
#'
#' Continuous covariates are assembled into a numeric matrix named `X`.
#' Categorical covariates are not included in this matrix; their integer
#' indices and level counts are prepared separately by
#' [create_insect_nimble_constants()]. When a temporal effect is requested,
#' the function validates and adds the temporal predictor. Coordinate-based
#' spatial effects use standardized longitude and latitude, whereas spatial
#' autocorrelation uses the original geographic coordinates in the modelling
#' data for subsequent processing by the constants-preparation function.
#'
#' The returned data list includes a `model_metadata` attribute containing
#' the response specification, covariate classifications, selected temporal
#' and spatial effects, coordinate column names, and coordinate
#' standardization parameters where applicable.
#'
#' @param modelling_data A data frame containing the observations used to
#'   construct the model data. It must contain the selected response columns,
#'   all specified covariate columns, and any required coordinate columns.
#'
#' @param modelCovs A character vector specifying the covariates included in
#'   the model. If omitted or `NULL`, no covariates are included. Duplicate
#'   covariate names are not allowed.
#'
#' @param covariate_types An optional named character vector identifying each
#'   covariate as `"continuous"` or `"categorical"`. Names must correspond to
#'   elements of `modelCovs`. If `NULL`, all specified covariates are treated
#'   as continuous.
#'
#' @param response A character string specifying which response to prepare.
#'   Supported values are `"biomass"`, `"diversity"`, and `"both"`. Defaults
#'   to `"biomass"`.
#'
#' @param biomass_col A character string naming the biomass response column
#'   in `modelling_data`. Defaults to `"sum_wet_weight"`. Required when
#'   `response` is `"biomass"` or `"both"`.
#'
#' @param diversity_col A character string naming the diversity response
#'   column in `modelling_data`. Required when `response` is `"diversity"` or
#'   `"both"`. Defaults to `NULL`.
#'
#' @param biomass_distribution A character string specifying the biomass
#'   observation distribution. Supported values are `"normal"` and
#'   `"lognormal"`. Defaults to `"normal"`. The log-normal distribution
#'   requires strictly positive biomass observations.
#'
#' @param diversity_distribution A character string specifying the diversity
#'   observation distribution. Supported values are `"poisson"` and
#'   `"negative_binomial"`. Defaults to `"poisson"`. Both distributions
#'   require non-negative integer-valued observations.
#'
#' @param temporal_data An optional numeric vector containing the temporal
#'   predictor for each observation. Required when `temporal_effect` is
#'   `"linear"` or `"cubic_cyclic"`. Its length must equal the number of rows
#'   in `modelling_data`, and its values must be finite and non-missing.
#'
#' @param temporal_effect A character string specifying the temporal
#'   structure. Supported values are `"none"`, `"linear"`, and
#'   `"cubic_cyclic"`. Defaults to `"none"`. The temporal basis and penalty
#'   matrix for a cyclic spline must be prepared separately and supplied to
#'   [create_insect_nimble_model()].
#'
#' @param spatial_effect A character string specifying the spatial structure.
#'   Supported values are `"none"`, `"coordinates"`, and `"autocorrelation"`.
#'   Defaults to `"none"`.
#'
#' @param longitude_col A character string naming the longitude column in
#'   `modelling_data`. Defaults to `"longitude"`. Required when
#'   `spatial_effect` is `"coordinates"` or `"autocorrelation"`.
#'
#' @param latitude_col A character string naming the latitude column in
#'   `modelling_data`. Defaults to `"latitude"`. Required when
#'   `spatial_effect` is `"coordinates"` or `"autocorrelation"`.
#'
#' @return A named list containing the data variables required by the
#'   selected model configuration. Depending on the arguments, the list
#'   contains:
#' \describe{
#'   \item{biomass}{Numeric vector of biomass observations, if biomass is
#'     included.}
#'   \item{diversity}{Numeric vector of diversity observations, if diversity
#'     is included.}
#'   \item{X}{Numeric design matrix containing the continuous covariates.
#'     Its rows correspond to observations and its columns correspond to
#'     `continuous_covariates`. If no continuous covariates are specified,
#'     an `N`-by-zero matrix is returned.}
#'   \item{temporal}{Numeric temporal predictor, if a temporal effect is
#'     requested.}
#'   \item{longitude, latitude}{Standardized longitude and latitude vectors
#'     when `spatial_effect = "coordinates"`. These are not added for
#'     `spatial_effect = "autocorrelation"`.}
#' }
#'
#' The returned list also has a `model_metadata` attribute containing:
#' \describe{
#'   \item{response}{Selected response specification.}
#'   \item{biomass_col, diversity_col}{Names of the selected response columns.}
#'   \item{biomass_distribution, diversity_distribution}{Selected response
#'     distributions.}
#'   \item{continuous_covariates}{Names of continuous covariates.}
#'   \item{categorical_covariates}{Names of categorical covariates.}
#'   \item{temporal_effect}{Selected temporal structure.}
#'   \item{spatial_effect}{Selected spatial structure.}
#'   \item{longitude_col, latitude_col}{Names of the coordinate columns.}
#'   \item{spatial_coordinate_center}{Named vector containing the longitude
#'     and latitude means used for standardization, or `NULL` when no
#'     coordinate-based effect is used.}
#'   \item{spatial_coordinate_scale}{Named vector containing the longitude
#'     and latitude standard deviations used for standardization, or `NULL`
#'     when no coordinate-based effect is used.}
#' }
#'
#' @details
#' The number of observations is determined from the rows of
#' `modelling_data`. Response variables and continuous covariates are checked
#' for numeric types, missing values, and non-finite values. The function
#' rejects non-positive biomass values for a log-normal model and rejects
#' negative or non-integer diversity observations for Poisson and
#' negative-binomial models.
#'
#' Categorical covariates are classified but are not converted to integer
#' indices in this function. Their indices and numbers of levels are created
#' by [create_insect_nimble_constants()].
#'
#' For `spatial_effect = "coordinates"`, longitude and latitude are
#' standardized independently by subtracting their sample means and dividing
#' by their sample standard deviations. Both coordinates must have
#' non-zero variation. The means and standard deviations are stored in the
#' `model_metadata` attribute so that the same transformation can be applied
#' consistently to new data.
#'
#' For `spatial_effect = "autocorrelation"`, coordinates are validated but
#' are not added as standardized predictors to the returned data list.
#' Spatial indices and the distance matrix are prepared separately by
#' [create_insect_nimble_constants()].
#'
#' For `temporal_effect = "cubic_cyclic"`, this function only prepares the
#' temporal predictor. Construction of the spline basis, penalty matrix,
#' and penalized eigenbasis is handled separately by
#' [create_insect_nimble_model()].
#'
#' @examples
#' \dontrun{
#' # Prepare data for a biomass model with continuous covariates.
#' nimble_data <- create_insect_nimble_data(
#'   modelling_data = insect_data,
#'   modelCovs = c("temperature", "precipitation"),
#'   response = "biomass",
#'   biomass_col = "sum_wet_weight",
#'   biomass_distribution = "lognormal"
#' )
#'
#' # Inspect the prepared response and design matrix.
#' str(nimble_data)
#' dim(nimble_data$X)
#'
#' # Prepare data for a joint biomass-diversity model.
#' joint_data <- create_insect_nimble_data(
#'   modelling_data = insect_data,
#'   modelCovs = c("temperature", "habitat"),
#'   covariate_types = c(
#'     temperature = "continuous",
#'     habitat = "categorical"
#'   ),
#'   response = "both",
#'   biomass_col = "sum_wet_weight",
#'   diversity_col = "species_richness",
#'   biomass_distribution = "lognormal",
#'   diversity_distribution = "negative_binomial"
#' )
#'
#' # Inspect metadata attached to the data list.
#' attr(joint_data, "model_metadata")
#' }
#'
#' @seealso
#' [create_insect_nimble_model()],
#' [create_insect_nimble_constants()]
#'
#' @export


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