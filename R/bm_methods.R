#' Methods for `bm_fit` objects
#'
#' `print()` gives a short description of the model and MCMC run. `summary()`
#' returns posterior summaries with R-hat and effective sample size for all
#' parameters except the individual random-effect levels (`u_*`), which can be
#' requested with `random_levels = TRUE`.
#'
#' @param x,object A `bm_fit` object from [fit_model()].
#' @param random_levels Include individual random-effect levels in `summary()`.
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `summary()` an object of class
#'   `summary.bm_fit` containing the table in `$table`.
#' @name bm_fit-methods
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 10, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", n_chains = 2, n_iter = 2000,
#'                  n_burnin = 1000, n_thin = 2, seed = 1, parallel = FALSE)
#' print(fit)
#' summary(fit)
#' }
#' @md
NULL

#' @rdname bm_fit-methods
#' @export
#' @md
print.bm_fit <- function(x, ...) {
  sp <- x$spec
  cat("<bm_fit> Bayesian hierarchical model\n")
  cat("  Responses :", paste0(sp$responses, " (", sp$family[sp$responses], ")",
                              collapse = ", "), "\n")
  fx <- vapply(sp$fixed_effects, function(v) {
    if (length(v)) paste(v, collapse = ", ") else "none"
  }, character(1))
  cat("  Fixed     :", paste0(names(fx), ": ", fx, collapse = "; "), "\n")
  cat("  Random    :", if (length(sp$random_effects)) {
    paste(sp$random_effects, collapse = ", ") } else "none", "\n")
  cat("  Temporal  :", sp$temporal_effect, "| Spatial:", sp$spatial_effect,
      if (isTRUE(sp$spatiotemporal_interaction)) "| with interaction", "\n")
  sh <- sp$sharing
  if (!is.null(sh) && bm_any_shared(sh)) {
    parts <- c(
      if (length(sh$fixed)) paste0("fixed (", paste(sh$fixed, collapse = ", "), ")"),
      if (length(sh$random)) paste0("random (", paste(sh$random, collapse = ", "), ")"),
      if (sh$spatial) "spatial", if (sh$temporal) "temporal")
    cat("  Shared    :", paste(parts, collapse = "; "),
        if (length(parts)) paste0("[", if (sh$loadings) "latent + loading" else "identical",
                                  "; reference: ", sh$ref, "]"), "\n")
    if (length(sh$correlated)) {
      cat("  Correlated:", paste(sh$correlated, collapse = ", "), "\n")
    }
  }
  m <- x$mcmc
  cat(sprintf("  MCMC      : %d chains x %d iter (burn-in %d, thin %d), %d saved/chain\n",
              m$n_chains, m$n_iter, m$n_burnin, m$n_thin,
              nrow(x$samples[[1]])))
  cat("  Chain seeds:", paste(x$seeds$chains, collapse = ", "), "\n")
  d <- x$diagnostics
  if (!all(is.na(d$rhat))) {
    cat(sprintf("  Convergence: max R-hat %.3f, min ESS %.0f\n",
                max(d$rhat, na.rm = TRUE), min(d$ess, na.rm = TRUE)))
  }
  invisible(x)
}

#' @rdname bm_fit-methods
#' @export
#' @md
summary.bm_fit <- function(object, random_levels = FALSE, ...) {
  d <- object$diagnostics
  if (!random_levels) d <- d[!grepl("^(u|w_cyc|w_sp|v_sp)_", d$parameter), ]
  structure(list(table = d, fit = object), class = "summary.bm_fit")
}

#' @rdname bm_fit-methods
#' @export
#' @md
print.summary.bm_fit <- function(x, digits = 3, ...) {
  print(x$fit)
  cat("\nPosterior summaries (parameters on the untruncated/model scale):\n")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = digits)
  print(tab, row.names = FALSE)
  invisible(x)
}
