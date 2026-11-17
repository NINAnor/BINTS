#' Residuals of a fitted biomass-richness model
#'
#' Observed-minus-fitted residuals at the observations, in three flavours. All
#' are conditional on the estimated random effects and spatial/temporal terms
#' (the same conditioning as `predict(fit)`), and are averaged over posterior
#' draws.
#'
#' * `"raw"`: observed minus the posterior mean of the **expected response**
#'   (the mean of the truncated distribution, see [predict.bm_fit()]), in the
#'   units of the response.
#' * `"standardised"`: Pearson residuals \eqn{(y - E)/\sqrt{V}} using the mean
#'   and variance of the truncated distribution in each posterior draw,
#'   averaged over draws. Approximately unit variance for a well specified
#'   model.
#' * `"posterior_predictive"`: randomised-quantile type residuals
#'   \eqn{\Phi^{-1}(u_i)} with \eqn{u_i = \overline{F(y_i \mid \theta)}}, the
#'   posterior predictive probability integral transform computed with the
#'   fitted (truncated) cdf. Approximately standard normal for a well specified
#'   model, and comparable across responses and families. For richness, which
#'   is a count, the continuous cdf is used.
#'
#' @param object A `bm_fit` object.
#' @param type `"raw"`, `"standardised"` or `"posterior_predictive"`.
#' @param responses Responses (default all fitted).
#' @param n_draws Number of posterior draws used.
#' @param random Random-effect mode; see [predict.bm_fit()].
#' @param seed Optional seed (the RNG state of the session is restored).
#' @param ... Unused.
#' @return A long data frame with `.row` (row of the original data), `response`,
#'   `observed`, `fitted` (posterior mean expected response) and `residual`.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 12, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", n_chains = 2, n_iter = 2000,
#'                  n_burnin = 1000, n_thin = 2, seed = 1, parallel = FALSE)
#' head(residuals(fit, type = "posterior_predictive", n_draws = 200))
#' }
#' @md
residuals.bm_fit <- function(object, type = c("raw", "standardised", "posterior_predictive"),
                             responses = NULL, n_draws = 500, random = "conditional",
                             seed = NULL, ...) {
  type <- match.arg(type)
  core <- bm_predict_core(object, NULL, responses, random, "simulate", "conditional",
                          n_draws, seed)
  fams <- bm_families()
  parts <- lapply(core$responses, function(r) {
    fam <- fams[[object$spec$family[[r]]]]
    mu <- core$mu[[r]]
    sg <- core$sigma[[r]]
    y <- object$prep$y[[r]]
    Y <- matrix(y, nrow(mu), length(y), byrow = TRUE)
    E <- fam$expected(mu, sg)
    res <- switch(
      type,
      raw = y - colMeans(E),
      standardised = colMeans((Y - E) / sqrt(fam$variance(mu, sg))),
      posterior_predictive = {
        u <- if (isTRUE(fam$discrete)) {
          colMeans((fam$cdf(Y, mu, sg) + fam$cdf_left(Y, mu, sg)) / 2)
        } else colMeans(fam$cdf(Y, mu, sg))
        stats::qnorm(pmin(pmax(u, 1e-10), 1 - 1e-10))
      }
    )
    data.frame(.row = object$prep$kept_rows, response = r, observed = y,
               fitted = colMeans(E), residual = res)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  attr(out, "type") <- type
  out
}
