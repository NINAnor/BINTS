run_insect_nimble_model <- function(



    modelling_data,



    modelCovs,



    covariate_types = NULL,





    # ==========================================================

    # Responses

    # ==========================================================



    response = "both",



    biomass_col = "sum_wet_weight",



    diversity_col = "shannon_div",





    # ==========================================================

    # Response distributions

    # ==========================================================



    biomass_distribution = "normal",



    diversity_distribution = "normal",





    # ==========================================================

    # Replicate structure

    # ==========================================================



    replicate_col = "sampling_name",





    # ==========================================================

    # Temporal structure

    # ==========================================================



    temporal_col = c(

        "start_date_julian",

        "end_date_julian"

    ),





    # ==========================================================

    # Random effects

    # ==========================================================



    random_year = TRUE,



    random_locality = TRUE,



    random_replicate = TRUE,





    # ==========================================================

    # Temporal effect

    # ==========================================================



    temporal_effect = "cubic_cyclic",


    temporal_basis_k = 10,





    # ==========================================================

    # Spatial effect

    #

    # "none"

    # "coordinates"

    # "autocorrelation"

    # ==========================================================



    spatial_effect = "none",



    # If TRUE:

    #   - coordinates: one shared spatial coordinate effect

    #   - autocorrelation: one shared latent spatial process

    #

    # If FALSE:

    #   - biomass and diversity receive separate spatial effects

    spatial_shared = TRUE,



    longitude_col = "longitude",



    latitude_col = "latitude",





    # ==========================================================

    # Spatial priors

    #

    # Range is specified in metres because coordinates are

    # transformed to EPSG:25833 for the distance calculation.

    # ==========================================================



    spatial_range_prior = c(

        1000,

        500000

    ),



    spatial_sigma_prior = 1,





    # ==========================================================

    # Output

    # ==========================================================



    output_dir = temp_dir,





    # ==========================================================

    # MCMC

    # ==========================================================



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



    if (!requireNamespace(

        "nimble",

        quietly = TRUE

    )) {



        stop(

            "Package 'nimble' is required.",

            call. = FALSE

        )



    }





    # ==========================================================

    # Validate response

    # ==========================================================



    response <- match.arg(

        response,

        choices = c(

            "biomass",

            "diversity",

            "both"

        )

    )





    # ==========================================================

    # Validate distributions

    # ==========================================================



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

            "normal",

            "poisson",

            "negative_binomial"

        )

    )





    # ==========================================================

    # Validate spatial effect

    # ==========================================================



    spatial_effect <- match.arg(

        spatial_effect,

        choices = c(

            "none",

            "coordinates",

            "autocorrelation"

        )

    )





    # ==========================================================

    # Validate spatial priors

    # ==========================================================



    if (

        length(spatial_range_prior) != 2 ||

        any(!is.finite(spatial_range_prior)) ||

        spatial_range_prior[1] <= 0 ||

        spatial_range_prior[2] <= spatial_range_prior[1]

    ) {



        stop(

            "'spatial_range_prior' must contain two finite ",

            "positive values with upper > lower.",

            call. = FALSE

        )



    }





    if (

        length(spatial_sigma_prior) != 1 ||

        !is.finite(spatial_sigma_prior) ||

        spatial_sigma_prior <= 0

    ) {



        stop(

            "'spatial_sigma_prior' must be a single positive ",

            "finite value.",

            call. = FALSE

        )



    }





    # ==========================================================

    # Validate modelling data

    # ==========================================================



    if (!is.data.frame(modelling_data)) {



        stop(

            "'modelling_data' must be a data.frame.",

            call. = FALSE

        )



    }





    if (nrow(modelling_data) == 0) {



        stop(

            "'modelling_data' contains no rows.",

            call. = FALSE

        )



    }





    # ==========================================================

    # Validate model covariates

    # ==========================================================



    missing_covariates <- setdiff(

        modelCovs,

        names(modelling_data)

    )



    if (length(missing_covariates) > 0) {



        stop(

            "The following model covariates were not found in ",

            "'modelling_data': ",

            paste(

                missing_covariates,

                collapse = ", "

            ),

            call. = FALSE

        )



    }





    # ==========================================================

    # Validate covariate types

    # ==========================================================



    if (is.null(covariate_types)) {



        covariate_types <- setNames(

            rep(

                "continuous",

                length(modelCovs)

            ),

            modelCovs

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





    invalid_types <- setdiff(

        unique(

            covariate_types[modelCovs]

        ),

        c(

            "continuous",

            "categorical"

        )

    )



    if (length(invalid_types) > 0) {



        stop(

            "Invalid covariate type(s): ",

            paste(

                invalid_types,

                collapse = ", "

            ),

            ". Allowed types are 'continuous' and 'categorical'.",

            call. = FALSE

        )



    }





    # ==========================================================

    # Validate spatial coordinates

    #

    # Only required when a spatial effect is requested.

    # ==========================================================



    if (spatial_effect != "none") {



        required_spatial_columns <- c(

            longitude_col,

            latitude_col

        )



        missing_spatial_columns <- setdiff(

            required_spatial_columns,

            names(modelling_data)

        )



        if (length(missing_spatial_columns) > 0) {



            stop(

                "The following spatial coordinate columns were ",

                "not found in 'modelling_data': ",

                paste(

                    missing_spatial_columns,

                    collapse = ", "

                ),

                call. = FALSE

            )



        }



        if (

            !is.numeric(modelling_data[[longitude_col]]) ||

            !is.numeric(modelling_data[[latitude_col]])

        ) {



            stop(

                "Spatial coordinate columns must be numeric.",

                call. = FALSE

            )



        }



        if (

            any(!is.finite(

                modelling_data[[longitude_col]]

            )) ||

            any(!is.finite(

                modelling_data[[latitude_col]]

            ))

        ) {



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





    # ----------------------------------------------------------

    # Biomass

    # ----------------------------------------------------------



    if (response %in% c(

        "biomass",

        "both"

    )) {



        if (!biomass_col %in% names(modelling_data)) {



            stop(

                "Biomass column not found: ",

                biomass_col,

                call. = FALSE

            )



        }



        biomass <- modelling_data[[biomass_col]]



    }





    # ----------------------------------------------------------

    # Diversity

    # ----------------------------------------------------------



    if (response %in% c(

        "diversity",

        "both"

    )) {



        if (is.null(diversity_col)) {



            stop(

                "diversity_col must be supplied when ",

                "response includes diversity.",

                call. = FALSE

            )



        }



        if (!diversity_col %in% names(modelling_data)) {



            stop(

                "Diversity column not found: ",

                diversity_col,

                call. = FALSE

            )



        }



        diversity <- modelling_data[[diversity_col]]



    }





    # ==========================================================

    # Replicate number

    # ==========================================================



    replicates <- NULL



    if (!is.null(replicate_col)) {



        if (!replicate_col %in% names(modelling_data)) {



            stop(

                "Replicate column not found: ",

                replicate_col,

                call. = FALSE

            )



        }



        required_replicate_cols <- c(

            "locality",

            "year"

        )



        missing_replicate_cols <- setdiff(

            required_replicate_cols,

            names(modelling_data)

        )



        if (length(missing_replicate_cols) > 0) {



            stop(

                "The following columns are required to construct ",

                "replicate IDs: ",

                paste(

                    missing_replicate_cols,

                    collapse = ", "

                ),

                call. = FALSE

            )



        }



        replicates <- modelling_data |>

            dplyr::mutate(



                replicate_number = as.numeric(

                    sub(

                        ".*-",

                        "",

                        .data[[replicate_col]]

                    )

                )



            ) |>

            dplyr::group_by(

                locality,

                year

            ) |>

            dplyr::mutate(



                sampling_id = as.numeric(

                    as.factor(

                        replicate_number

                    )

                )



            ) |>

            dplyr::ungroup() |>

            dplyr::pull(

                sampling_id

            )



    }



# ==========================================================
# Temporal effect
# ==========================================================

temporal_data <- NULL

# These are retained for the NIMBLE constants if needed
temporal_basis_penalized <- NULL
temporal_penalty_eigenvalues <- NULL
n_temporal_basis <- 0L

# Original mgcv objects
temporal_basis <- NULL
temporal_penalty <- NULL

temporal_effect <- match.arg(
    temporal_effect,
    choices = c(
        "none",
        "linear",
        "cubic_cyclic"
    )
)

if (!is.numeric(temporal_basis_k) ||
    length(temporal_basis_k) != 1L ||
    !is.finite(temporal_basis_k) ||
    temporal_basis_k < 4 ||
    temporal_basis_k != as.integer(temporal_basis_k)) {

    stop(
        "'temporal_basis_k' must be a single integer >= 4.",
        call. = FALSE
    )
}

temporal_basis_k <- as.integer(
    temporal_basis_k
)

if (temporal_effect != "none") {

    if (is.null(temporal_col)) {

        stop(
            "temporal_col must be supplied when ",
            "temporal_effect is not 'none'.",
            call. = FALSE
        )

    }

    missing_temporal <- setdiff(
        temporal_col,
        names(modelling_data)
    )

    if (length(missing_temporal) > 0) {

        stop(
            "Temporal columns not found: ",
            paste(
                missing_temporal,
                collapse = ", "
            ),
            call. = FALSE
        )

    }

    temporal_data <- rowMeans(
        modelling_data[
            ,
            temporal_col,
            drop = FALSE
        ],
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

        temporal_data <- as.numeric(
            scale(temporal_data)
        )

    }

    # ------------------------------------------------------
    # Cubic cyclic temporal effect
    # ------------------------------------------------------

    if (temporal_effect == "cubic_cyclic") {

        if (!requireNamespace(
            "mgcv",
            quietly = TRUE
        )) {

            stop(
                "Package 'mgcv' is required for ",
                "temporal_effect = 'cubic_cyclic'.",
                call. = FALSE
            )

        }

        # IMPORTANT:
        # Do NOT standardise temporal_data here.
        #
        # For a cyclic spline, the original temporal scale
        # defines the cyclic boundary.
        temporal_smooth_data <- data.frame(
            temporal = as.numeric(
                temporal_data
            )
        )

        temporal_smooth_spec <- mgcv::s(
            temporal,
            bs = "cc",
            k = temporal_basis_k
        )

        # Explicit cyclic boundary.
        temporal_knots <- list(
            temporal = range(
                temporal_smooth_data$temporal,
                finite = TRUE
            )
        )

        temporal_smooth <- mgcv::smoothCon(
            temporal_smooth_spec,
            data = temporal_smooth_data,
            knots = temporal_knots,
            absorb.cons = FALSE,
            scale.penalty = FALSE
        )[[1]]

        # --------------------------------------------------
        # Save the ORIGINAL mgcv basis and penalty.
        #
        # create_insect_nimble_model() will perform the
        # eigen-decomposition and construct the penalised
        # eigenbasis.
        # --------------------------------------------------

        temporal_basis <- temporal_smooth$X

        temporal_penalty <- temporal_smooth$S[[1]]

        if (!is.matrix(temporal_basis) ||
            !is.numeric(temporal_basis)) {

            stop(
                "Failed to construct a valid mgcv temporal basis.",
                call. = FALSE
            )

        }

        if (!is.matrix(temporal_penalty) ||
            !is.numeric(temporal_penalty)) {

            stop(
                "Failed to construct a valid mgcv temporal penalty.",
                call. = FALSE
            )

        }

        if (any(!is.finite(temporal_basis)) ||
            any(!is.finite(temporal_penalty))) {

            stop(
                "The temporal basis or penalty contains ",
                "non-finite values.",
                call. = FALSE
            )

        }

    }

}






    # ==========================================================

    # Identify continuous covariates

    # ==========================================================



    continuous_covariates <- modelCovs[

        covariate_types[modelCovs] == "continuous"

    ]





    # ==========================================================

    # Scale continuous covariates

    # ==========================================================



    if (length(continuous_covariates) > 0) {



        scale_covariate <- function(x) {



            if (!is.numeric(x)) {



                stop(

                    "Continuous covariates must be numeric.",

                    call. = FALSE

                )



            }



            mean_x <- mean(

                x,

                na.rm = TRUE

            )



            sd_x <- sd(

                x,

                na.rm = TRUE

            )



            if (

                !is.finite(sd_x) ||

                sd_x == 0

            ) {



                stop(

                    "A continuous covariate has zero or undefined ",

                    "standard deviation.",

                    call. = FALSE

                )



            }



            as.numeric(

                (x - mean_x) / sd_x

            )



        }





        modelling_data <- modelling_data |>

            dplyr::mutate(

                dplyr::across(

                    dplyr::all_of(

                        continuous_covariates

                    ),

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

    # Supply the original mgcv basis and penalty
    temporal_basis = temporal_basis,

    temporal_penalty = temporal_penalty,

    output_dir = output_dir,

    covariate_types = covariate_types,

    response = response,

    biomass_distribution = biomass_distribution,

    diversity_distribution = diversity_distribution,

    # ------------------------------------------------------
    # Spatial structure
    # ------------------------------------------------------

    spatial_effect = spatial_effect,

    spatial_shared = spatial_shared

)




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

    # ------------------------------------------------------
    # Spatial structure
    # ------------------------------------------------------

    spatial_effect = spatial_effect,

    longitude_col = longitude_col,

    latitude_col = latitude_col

)




    # ==========================================================

    # Create NIMBLE constants

    #

    # Important:

    # For autocorrelation, this function calculates ONLY the

    # fixed distance matrix. The spatial covariance matrix is

    # constructed inside NIMBLE from estimated rho and sigma.

    # ==========================================================



    constants <- create_insect_nimble_constants(



        modelling_data = modelling_data,



        modelCovs = modelCovs,



        covariate_types = covariate_types,



        replicate_col = replicate_col,



        random_year = random_year,



        random_locality = random_locality,



        random_replicate = random_replicate,



        # ------------------------------------------------------

        # Spatial structure

        # ------------------------------------------------------



        spatial_effect = spatial_effect,



        longitude_col = longitude_col,



        latitude_col = latitude_col



    )





    # ==========================================================

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


    # Create initial values

    # ==========================================================



    initial_values <- create_insect_nimble_initial_values(



        modelling_data = modelling_data,



        modelCovs = modelCovs,



        covariate_types = covariate_types,



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



        seed = seed,



        # ------------------------------------------------------

        # Spatial structure

        # ------------------------------------------------------



        spatial_effect = spatial_effect,



        spatial_shared = spatial_shared,



        longitude_col = longitude_col,



        latitude_col = latitude_col,



        spatial_range_prior = spatial_range_prior,



        spatial_sigma_prior = spatial_sigma_prior



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





    source(

        model_file$file,

        local = TRUE

    )





    model_code <- biomass_diversity_model





    # ==========================================================

    # Build NIMBLE model

    # ==========================================================



    nimble_model <- nimble::nimbleModel(



        code = model_code,



        constants = constants,



        data = data,



        inits = initial_values,



        check = TRUE,



        calculate = FALSE



    )



    nimble_model$initializeInfo()



    # ==========================================================

    # Compile NIMBLE model

    # ==========================================================



    compiled_model <- nimble::compileNimble(

        nimble_model

    )





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
                spatial_shared = spatial_shared
            )
        )
    }


    # ==========================================================
    # Configure MCMC
    # ==========================================================

    mcmc_conf <- nimble::configureMCMC(
        nimble_model,
        enableWAIC = run_waic
    )

    # ==========================================================
    # Build MCMC
    # ==========================================================

    mcmc <- nimble::buildMCMC(
        mcmc_conf
    )

    # ==========================================================
    # Compile MCMC
    # ==========================================================

    compiled_mcmc <- nimble::compileNimble(
        mcmc,
        project = nimble_model
    )

    # ==========================================================
    # Run MCMC
    # ==========================================================
    samples <- nimble::runMCMC(

        compiled_mcmc,

        niter = niter,

        nburnin = nburnin,

        thin = thin,

        nchains = nchains,

        setSeed = TRUE,

        samplesAsCodaMCMC = TRUE,

        summary = TRUE,

        WAIC = run_waic
    )


    # ==========================================================

    # Return everything

    # ==========================================================

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
        temporal_basis_k = temporal_basis_k,
        temporal_basis_penalized = temporal_basis_penalized,
        temporal_penalty_eigenvalues = temporal_penalty_eigenvalues,

        modelCovs = modelCovs,

        covariate_types = covariate_types,

        # ------------------------------------------------------

        # Spatial model information

        # ------------------------------------------------------

        spatial_effect = spatial_effect,

        spatial_shared = spatial_shared,

        longitude_col = longitude_col,

        latitude_col = latitude_col,

        spatial_range_prior = spatial_range_prior,

        spatial_sigma_prior = spatial_sigma_prior

    )

}