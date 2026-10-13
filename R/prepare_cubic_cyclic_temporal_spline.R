
# ============================================================================
# Prepare cubic cyclic temporal spline basis
# ============================================================================
#
# Constructs the mgcv cyclic cubic spline basis and penalty matrix that are
# supplied to create_insect_nimble_model().
#
# The returned object contains:
#
#   basis   : observation-by-basis matrix
#   penalty : mgcv spline penalty matrix
#   knots   : cyclic boundary knots
#
# Example:
#
#   temporal_spline <- prepare_cubic_cyclic_temporal_spline(
#     temporal = dat$temporal,
#     k = 10,
#     period = 1
#   )
#
#   model <- create_insect_nimble_model(
#     modelCovs = modelCovs,
#     temporal_effect = "cubic_cyclic",
#     temporal_basis = temporal_spline$basis,
#     temporal_penalty = temporal_spline$penalty
#   )
#
# For a cyclic effect, the temporal variable should represent a quantity
# whose endpoints are equivalent, e.g. day-of-year or another periodic phase.
# ============================================================================

prepare_cubic_cyclic_temporal_spline <- function(

    temporal,

    k = 10,

    period = NULL

) {

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required to construct the cyclic cubic spline.",
      call. = FALSE
    )
  }

  if (!is.numeric(temporal) || any(!is.finite(temporal))) {
    stop(
      "'temporal' must be a numeric vector containing only finite values.",
      call. = FALSE
    )
  }

  if (length(temporal) < 2) {
    stop(
      "'temporal' must contain at least two observations.",
      call. = FALSE
    )
  }

  if (!is.numeric(k) || length(k) != 1 || k < 4 || k != as.integer(k)) {
    stop(
      "'k' must be a single integer >= 4.",
      call. = FALSE
    )
  }

  lower <- min(temporal)

  if (is.null(period)) {
    upper <- max(temporal)
  } else {
    if (!is.numeric(period) || length(period) != 1 ||
        !is.finite(period) || period <= 0) {
      stop(
        "'period' must be a single positive finite value.",
        call. = FALSE
      )
    }
    upper <- lower + period
  }

  if (upper <= lower) {
    stop(
      "The temporal range must be positive.",
      call. = FALSE
    )
  }

  temporal_data <- data.frame(
    temporal = temporal
  )

  smooth <- mgcv::smoothCon(
    mgcv::s(
      temporal,
      bs = "cc",
      k = as.integer(k)
    ),
    data = temporal_data,
    knots = list(
      temporal = c(lower, upper)
    )
  )[[1]]

  list(
    basis = smooth$X,
    penalty = smooth$S[[1]],
    knots = c(lower, upper),
    k = ncol(smooth$X)
  )
}
