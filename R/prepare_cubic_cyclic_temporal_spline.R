#' Prepare a Cubic Cyclic Temporal Spline Basis
#'
#' Construct a cubic cyclic regression spline basis and its associated penalty
#' matrix using \pkg{mgcv}. The resulting basis and penalty can be supplied to
#' [create_insect_nimble_model()] to define a penalised cyclic temporal effect
#' in the NIMBLE insect biomass and diversity models.
#'
#' Cyclic splines are appropriate for periodic temporal variables whose
#' endpoints represent equivalent states, such as day of year. The function
#' defines the cyclic boundary knots from the observed temporal range or from
#' an explicitly specified period.
#'
#' @param temporal Numeric vector containing the temporal values for all
#'   observations. Must contain at least two finite values.
#' @param k Integer specifying the requested basis dimension for the cyclic
#'   cubic regression spline. Must be a single integer greater than or equal
#'   to 4. Defaults to `10`.
#' @param period Optional positive finite numeric value specifying the period
#'   of the temporal variable. If `NULL`, the upper boundary knot is the
#'   maximum observed temporal value. Otherwise, the upper boundary knot is
#'   calculated as `min(temporal) + period`.
#'
#' @details
#' The function uses [mgcv::smoothCon()] to construct a cyclic cubic regression
#' spline with `bs = "cc"`. The boundary knots are set explicitly using the
#' lower and upper temporal boundaries.
#'
#' When `period = NULL`, the boundary knots are the minimum and maximum values
#' of `temporal`. When a period is supplied, the lower boundary is the minimum
#' observed temporal value and the upper boundary is the lower boundary plus
#' the specified period.
#'
#' The returned basis matrix has one row per observation. The penalty matrix
#' is extracted from the first spline penalty returned by `mgcv`. The actual
#' number of basis columns is determined by the constructed spline and is
#' returned as `k`.
#'
#' The returned basis and penalty are intended for use with
#' [create_insect_nimble_model()]. Ensure that the temporal variable and
#' boundary knots correctly represent the intended cycle. In particular, when
#' supplying `period`, the observations should be compatible with the specified
#' cyclic domain.
#'
#' This function requires the \pkg{mgcv} package. It stops with an error if
#' `temporal` contains non-finite values, has fewer than two observations, or
#' if `k` or `period` fails the stated validation requirements.
#'
#' @return A named list with four components:
#' \describe{
#'   \item{`basis`}{Numeric matrix containing the constructed spline basis,
#'   with one row per observation.}
#'   \item{`penalty`}{Numeric matrix containing the spline penalty matrix
#'   returned by `mgcv`.}
#'   \item{`knots`}{Numeric vector of length two containing the lower and
#'   upper cyclic boundary knots.}
#'   \item{`k`}{Integer giving the number of columns in the constructed basis
#'   matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' temporal_spline <- prepare_cubic_cyclic_temporal_spline(
#'   temporal = dat$temporal,
#'   k = 10,
#'   period = 1
#' )
#'
#' model <- create_insect_nimble_model(
#'   modelCovs = modelCovs,
#'   temporal_effect = "cubic_cyclic",
#'   temporal_basis = temporal_spline$basis,
#'   temporal_penalty = temporal_spline$penalty
#' )
#' }
#'
#' @seealso
#' [mgcv::smoothCon()], [mgcv::s()],
#' [create_insect_nimble_model()]
#'
#' @importFrom mgcv s smoothCon
#' @export

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
