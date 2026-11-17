#' Observed versus predicted values
#'
#' Scatter plot of observed (x) against predicted (y) values at the
#' observations with the 1:1 line, optional credible (or prediction) intervals,
#' and RMSE and MAE per panel. Predictions are the posterior mean (or median)
#' of the **expected response** (the truncated mean, see [predict.bm_fit()])
#' conditional on the estimated random effects, i.e. fitted values.
#'
#' With `type = "predictive"` the predicted value is the centre of the
#' posterior predictive distribution and the intervals are prediction
#' intervals (they include observation noise, and are wider). The metrics
#' always use the plotted point prediction: `observed - predicted`.
#'
#' Panels have equal x and y limits so that the 1:1 line is the diagonal.
#' To avoid repeating predictions, compute `predict(model)` once and pass it as
#' `predictions` (the result for `newdata = NULL`).
#'
#' @param model A `bm_fit` object.
#' @param response One or more fitted responses; several are faceted.
#' @param type `"expected"` (default) or `"predictive"`.
#' @param statistic Point prediction: `"mean"` or `"median"` of the posterior.
#' @param interval Draw posterior intervals as vertical lines.
#' @param level Interval level.
#' @param facet_by `NULL`, `"response"`, or the name of a column of the modelling
#'   data such as `"year"` or `"locality"`.
#' @param metrics Annotate each panel with RMSE and MAE.
#' @param predictions Optional precomputed `predict(model, ...)` result for the
#'   original observations (long format with `.row` and `response`).
#' @param n_draws,random,seed Passed to [predict.bm_fit()] when predictions are
#'   computed here.
#' @param max_facets Upper limit on the number of panels.
#' @return A `ggplot2` object; its `$data` holds the plotted values.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 12, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", fixed_effects = "x1", n_chains = 2,
#'                  n_iter = 2000, n_burnin = 1000, n_thin = 2, seed = 1,
#'                  parallel = FALSE)
#' plot_observed_predicted(fit, "biomass", interval = TRUE, facet_by = "year")
#' }
#' @md
plot_observed_predicted <- function(model, response = "biomass",
                                    type = c("expected", "predictive"),
                                    statistic = c("mean", "median"), interval = FALSE,
                                    level = 0.95, facet_by = NULL, metrics = TRUE,
                                    predictions = NULL, n_draws = 500,
                                    random = "conditional", max_facets = 30, seed = NULL) {
  type <- match.arg(type)
  statistic <- match.arg(statistic)
  response <- bm_check_responses(model, response)
  if (!is.null(facet_by) && !(length(facet_by) == 1L &&
                              (identical(facet_by, "response") ||
                               facet_by %in% names(model$prep$data)))) {
    stop("`facet_by` must be NULL, 'response' or a column of the modelling data.",
         call. = FALSE)
  }
  pred <- if (is.null(predictions)) {
    stats::predict(model, responses = response, type = type, level = level,
                   n_draws = n_draws, random = random, seed = seed)
  } else predictions
  need <- c(".row", "response", "mean", "median", "lower", "upper")
  if (!is.data.frame(pred) || !all(need %in% names(pred))) {
    stop("`predictions` must be the output of predict() for the observations ",
         "(columns: ", paste(need, collapse = ", "), ").", call. = FALSE)
  }
  if (inherits(pred, "sf")) pred <- sf::st_drop_geometry(pred)
  pred <- as.data.frame(pred)
  pred <- pred[pred$response %in% response, , drop = FALSE]
  if (!nrow(pred)) stop("`predictions` contain none of the requested responses.", call. = FALSE)
  if (identical(attr(predictions, "observations"), FALSE) ||
      !all(pred$.row %in% model$prep$kept_rows)) {
    stop("`predictions` must be for the observations used in the fit (newdata = NULL).",
         call. = FALSE)
  }
  pred$observed <- NA_real_
  for (r in response) {
    i <- pred$response == r
    pred$observed[i] <- bm_obs_for(model, r, pred$.row[i])
  }
  pred$predicted <- pred[[statistic]]
  pred <- pred[stats::complete.cases(pred[c("observed", "predicted")]), , drop = FALSE]

  by_col <- !is.null(facet_by) && facet_by != "response"
  pred$.facet <- if (by_col) bm_data_col(model, facet_by, pred$.row) else "all"
  pred$.facet <- factor(pred$.facet)
  if (by_col) bm_check_facets(nlevels(pred$.facet) * length(response), max_facets, facet_by)

  ycols <- if (interval) c("predicted", "lower", "upper") else "predicted"
  blank <- bm_square_blank(pred, c("response", ".facet"), "observed", ycols)
  p <- ggplot2::ggplot(pred, ggplot2::aes(x = .data$observed, y = .data$predicted)) +
    ggplot2::geom_blank(data = blank, ggplot2::aes(x = .data$.lim, y = .data$.lim),
                        inherit.aes = FALSE) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "firebrick", linetype = 2)
  if (interval) {
    p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                     alpha = 0.25, colour = "grey40")
  }
  p <- p + ggplot2::geom_point(alpha = 0.5, size = 1.2)
  if (metrics) {
    key <- interaction(pred$response, pred$.facet, drop = TRUE)
    mt <- do.call(rbind, lapply(split(seq_len(nrow(pred)), key), function(ii) {
      e <- pred$observed[ii] - pred$predicted[ii]
      data.frame(response = pred$response[ii[1]], .facet = pred$.facet[ii[1]],
                 label = sprintf("RMSE = %s\nMAE = %s", format(sqrt(mean(e^2)), digits = 3),
                                 format(mean(abs(e)), digits = 3)))
    }))
    p <- p + ggplot2::geom_text(data = mt, ggplot2::aes(label = .data$label), x = -Inf,
                                y = Inf, hjust = -0.08, vjust = 1.3, size = 3,
                                inherit.aes = FALSE)
  }
  if (by_col && length(response) > 1L) {
    p <- p + ggplot2::facet_grid(response ~ .facet, scales = "free")
  } else if (by_col) {
    p <- p + ggplot2::facet_wrap(~ .facet, scales = "free")
  } else if (length(response) > 1L) {
    p <- p + ggplot2::facet_wrap(~ response, scales = "free")
  }
  ylab <- paste0("Predicted (", statistic, if (type == "expected") " expected response" else
    " posterior predictive", ")")
  p + ggplot2::labs(x = "Observed", y = ylab)
}

#' Residuals against covariates
#'
#' Plots residuals (see [residuals.bm_fit()]) against selected covariates, with
#' a horizontal reference line at zero and an optional smooth trend, in a grid
#' of panels by response and covariate. Residuals are observed minus predicted
#' (posterior mean expected response) for `type = "raw"`; standardised and
#' posterior predictive residuals are also available.
#'
#' Covariates can be any column of the modelling data, whether or not it is in the
#' model (for example `"longitude"` or `"latitude"`). Character, factor and
#' logical columns, and numeric columns with at most `max_numeric_levels`
#' distinct values (such as year), are treated as categorical and shown as
#' boxplots with jittered points; categorical covariates with more than
#' `max_levels` levels (such as locality) show the mean residual and an
#' approximate 95% interval per level. Continuous covariates get points and, if
#' `smooth = TRUE`, a penalised spline smooth.
#'
#' A single plot cannot mix continuous and discrete axes. When both kinds of
#' covariate are requested, categorical covariates are drawn at integer
#' positions in level order and their levels (or their number, if more than
#' four) are listed in the panel label.
#'
#' @param model A `bm_fit` object.
#' @param covariates Names of columns of the modelling data.
#' @param response One or more fitted responses.
#' @param type `"raw"`, `"standardised"` or `"posterior_predictive"`.
#' @param smooth Add a smooth trend for continuous covariates.
#' @param n_draws,random,seed Passed to [residuals.bm_fit()].
#' @param max_levels Categorical covariates with more levels use the summary
#'   display.
#' @param max_numeric_levels Numeric columns with at most this many distinct
#'   values are treated as categorical.
#' @return A `ggplot2` object.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 12, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", n_chains = 2, n_iter = 2000,
#'                  n_burnin = 1000, n_thin = 2, seed = 1, parallel = FALSE)
#' plot_residuals(fit, c("longitude", "latitude", "x1"), "biomass")
#' }
#' @md
plot_residuals <- function(model, covariates, response = "biomass",
                           type = c("raw", "standardised", "posterior_predictive"),
                           smooth = TRUE, n_draws = 500, random = "conditional",
                           max_levels = 30, max_numeric_levels = 10, seed = NULL) {
  type <- match.arg(type)
  response <- bm_check_responses(model, response)
  cinfo <- bm_cov_info(model$prep$data, covariates, max_numeric_levels)
  res <- stats::residuals(model, type = type, responses = response, n_draws = n_draws,
                          random = random, seed = seed)
  for (cv in covariates) res[[cv]] <- bm_data_col(model, cv, res$.row)
  res$.y <- res$residual
  long <- bm_cov_long(res, cinfo, c("response", ".y"))
  long <- long[stats::complete.cases(long[c(".x", ".y")]), , drop = FALSE]
  ylab <- switch(type, raw = "Residual (observed - predicted)",
                 standardised = "Standardised residual",
                 posterior_predictive = "Posterior predictive residual (normal scores)")
  p <- ggplot2::ggplot(long) +
    ggplot2::geom_hline(yintercept = 0, colour = "firebrick", linetype = 2) +
    bm_cov_layers(long, smooth, TRUE, max_levels)
  p <- p + if (length(response) > 1L) {
    ggplot2::facet_grid(response ~ .covariate, scales = "free",
                         labeller = ggplot2::label_wrap_gen(22))
  } else {
    ggplot2::facet_wrap(~ .covariate, scales = "free_x",
                         labeller = ggplot2::label_wrap_gen(22))
  }
  if (cinfo$mode == "factor" && any(vapply(cinfo$covariates, function(ci) {
    length(ci$levels) > max_levels }, logical(1)))) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  }
  p + ggplot2::labs(x = NULL, y = ylab)
}

# Average covariate effect curve for response r. The focal covariate is set to
# each grid value for a subsample of the observed rows, expected responses are
# predicted, and averaged over the rows within every posterior draw. All other
# covariates keep their observed values. "conditional" keeps the estimated
# random effects, spatial and temporal terms of each row; "marginal" integrates
# over the random effects (population level).
bm_covariate_effect <- function(model, r, cv, mode, ci, level, n_draws, n_rows, n_grid,
                                seed) {
  prep <- model$prep
  dat <- prep$data
  allowed <- c(model$spec$fixed_effects[[r]],
               if (identical(model$spec$temporal_effect, "cyclic")) model$spec$columns$season)
  if (!cv %in% allowed) {
    stop("Effect curves are available for fixed-effect covariates of '", r, "'",
         if (length(allowed)) paste0(" (", paste(allowed, collapse = ", "), ")"),
         " and the seasonal variable; '", cv, "' is not one.", call. = FALSE)
  }
  rows <- if (nrow(dat) > n_rows) unique(round(seq(1, nrow(dat), length.out = n_rows))) else
    seq_len(nrow(dat))
  base <- dat[rows, , drop = FALSE]
  grid <- if (ci$cat) ci$levels else {
    q <- stats::quantile(dat[[cv]], c(0.01, 0.99), names = FALSE)
    seq(q[1], q[2], length.out = n_grid)
  }
  nd <- do.call(rbind, lapply(seq_along(grid), function(k) {
    b <- base
    b[[cv]] <- if (is.numeric(dat[[cv]]) && ci$cat) as.numeric(grid[k]) else
      if (is.factor(dat[[cv]])) factor(grid[k], levels = levels(dat[[cv]])) else grid[k]
    b
  }))
  dr <- stats::predict(model, nd, responses = r, type = "expected", summarise = FALSE,
                       random = if (mode == "marginal") "marginal" else "conditional",
                       n_draws = n_draws, seed = seed)
  M <- dr$draws[[r]]
  blk <- rep(seq_along(grid), each = length(rows))
  avg <- vapply(seq_along(grid), function(k) rowMeans(M[, blk == k, drop = FALSE]),
                numeric(nrow(M)))
  a <- (1 - level) / 2
  out <- data.frame(
    response = r, covariate = cv, grid = if (ci$cat) grid else as.numeric(grid),
    median = apply(avg, 2, stats::median),
    lower = apply(avg, 2, stats::quantile, probs = a, names = FALSE),
    upper = apply(avg, 2, stats::quantile, probs = 1 - a, names = FALSE)
  )
  out
}

#' Observed and predicted responses against covariates
#'
#' For each response and covariate, plots the observed values and the fitted
#' (predicted) values in separate panels with the same y-axis limits, so they
#' can be compared directly, and a shared x-axis per covariate.
#'
#' @section What is shown:
#' * **Fitted predictions at observed covariate values** (always): the
#'   "Predicted" panel shows, for every observation, the posterior mean expected
#'   response from the fitted model, i.e. all terms (the other covariates, random
#'   effects, space and time) at their observed values. The spread against a
#'   single covariate therefore reflects all other terms, and it is not an
#'   effect of that covariate. `trend = TRUE` adds a penalised spline smooth
#'   (with its approximate 95% band) through the plotted points of each panel; it
#'   describes the points, it is not a model-based effect.
#' * **Covariate effect curves** (`effect`, optional): a model-based
#'   effect of one covariate, in the "Predicted" panel with a credible band. The
#'   focal covariate is set to each value of a grid for a subsample of rows and
#'   the predicted expected responses are averaged over those rows within every
#'   posterior draw (an average partial effect), keeping the other covariates at
#'   their observed values. `"conditional"` keeps the estimated random effects,
#'   spatial and temporal terms of each row. `"marginal"` integrates over the
#'   random effects instead (population level). Curves are available for
#'   fixed-effect covariates of the response (and the seasonal variable of a
#'   cyclic model); because the effect is on the response scale of a truncated
#'   distribution, it is not exactly linear in the coefficient.
#'
#' Continuous and categorical covariates are handled as in [plot_residuals()].
#'
#' @param model A `bm_fit` object.
#' @param covariates Names of columns of the modelling data.
#' @param response One or more fitted responses.
#' @param trend Add smooth trends through the plotted points.
#' @param effect `"none"` (default), `"conditional"` or `"marginal"`; see above.
#' @param level Credible level of effect curves.
#' @param n_draws Posterior draws used for fitted values and effect curves.
#' @param n_rows,n_grid Rows subsampled and grid values for effect curves.
#' @param random,seed Passed to prediction.
#' @param max_levels,max_numeric_levels As in [plot_residuals()].
#' @return A `ggplot2` object.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bm_data(n_year = 4, n_locality = 12, n_rep = 3, seed = 1)
#' fit <- fit_model(sim, responses = "biomass", fixed_effects = c("x1", "x2"),
#'                  n_chains = 2, n_iter = 2000, n_burnin = 1000, n_thin = 2,
#'                  seed = 1, parallel = FALSE)
#' plot_response_covariates(fit, c("x1", "year"), "biomass", effect = "marginal")
#' }
#' @md
plot_response_covariates <- function(model, covariates, response = "biomass", trend = TRUE,
                                     effect = c("none", "conditional", "marginal"),
                                     level = 0.95, n_draws = 300, n_rows = 300,
                                     n_grid = 25, random = "conditional",
                                     max_levels = 30, max_numeric_levels = 10,
                                     seed = NULL) {
  effect <- match.arg(effect)
  response <- bm_check_responses(model, response)
  cinfo <- bm_cov_info(model$prep$data, covariates, max_numeric_levels)
  pr <- stats::predict(model, responses = response, type = "expected", n_draws = n_draws,
                       random = random, seed = seed)
  pr <- as.data.frame(pr)
  pr$.observed <- NA_real_
  for (r in response) {
    i <- pr$response == r
    pr$.observed[i] <- bm_obs_for(model, r, pr$.row[i])
  }
  for (cv in covariates) pr[[cv]] <- bm_data_col(model, cv, pr$.row)
  src_obs <- pr
  src_obs$.y <- pr$.observed
  src_obs$.source <- "Observed"
  src_prd <- pr
  src_prd$.y <- pr$mean
  src_prd$.source <- "Predicted"
  src <- rbind(src_obs, src_prd)
  src$.source <- factor(src$.source, levels = c("Observed", "Predicted"))
  long <- bm_cov_long(src, cinfo, c("response", ".y", ".source"))
  long <- long[stats::complete.cases(long[c(".x", ".y")]), , drop = FALSE]

  eff <- NULL
  if (effect != "none") {
    eff <- do.call(rbind, lapply(response, function(r) {
      do.call(rbind, lapply(covariates, function(cv) {
        ci <- cinfo$covariates[[cv]]
        e <- bm_covariate_effect(model, r, cv, effect, ci, level, n_draws, n_rows, n_grid, seed)
        e$.x <- bm_x_value(e$grid, ci, cinfo$mode)
        e$.covariate <- bm_cov_label(ci, cinfo$mode)
        e$.cat <- ci$cat
        e$.source <- "Predicted"
        e
      }))
    }))
    eff$.covariate <- factor(eff$.covariate, levels = levels(long$.covariate))
    eff$.source <- factor(eff$.source, levels = levels(long$.source))
    if (cinfo$mode == "factor") eff$.x <- factor(as.character(eff$.x), levels = levels(long$.x))
  }

  # Same y limits for observed and predicted panels of a response
  lim <- do.call(rbind, lapply(response, function(r) {
    v <- c(long$.y[long$response == r], eff$lower[eff$response == r], eff$upper[eff$response == r])
    data.frame(response = r, .lim = range(v, na.rm = TRUE))
  }))
  first_x <- do.call(rbind, lapply(levels(long$.covariate), function(cl) {
    x <- long$.x[long$.covariate == cl][1]
    data.frame(.covariate = cl, .x = I(list(x)))
  }))
  blank <- do.call(rbind, lapply(seq_len(nrow(lim)), function(k) {
    do.call(rbind, lapply(levels(long$.covariate), function(cl) {
      do.call(rbind, lapply(levels(long$.source), function(s) {
        data.frame(response = lim$response[k], .covariate = cl, .source = s,
                   .y = lim$.lim[k], .x = first_x$.x[[match(cl, first_x$.covariate)]])
      }))
    }))
  }))
  blank$.covariate <- factor(blank$.covariate, levels = levels(long$.covariate))
  blank$.source <- factor(blank$.source, levels = levels(long$.source))
  if (cinfo$mode == "factor") blank$.x <- factor(as.character(blank$.x), levels = levels(long$.x))

  p <- ggplot2::ggplot(long) +
    ggplot2::geom_blank(data = blank, ggplot2::aes(x = .data$.x, y = .data$.y)) +
    bm_cov_layers(long, trend, TRUE, max_levels)
  if (!is.null(eff)) {
    ec <- eff[eff$.cat, , drop = FALSE]
    en <- eff[!eff$.cat, , drop = FALSE]
    if (nrow(en)) {
      p <- p + ggplot2::geom_ribbon(data = en, ggplot2::aes(x = .data$.x, ymin = .data$lower,
                                                            ymax = .data$upper),
                                    inherit.aes = FALSE, fill = "darkorange", alpha = 0.3) +
        ggplot2::geom_line(data = en, ggplot2::aes(x = .data$.x, y = .data$median),
                           colour = "darkorange", linewidth = 1)
    }
    if (nrow(ec)) {
      p <- p + ggplot2::geom_pointrange(data = ec, ggplot2::aes(x = .data$.x, y = .data$median,
                                                                ymin = .data$lower,
                                                                ymax = .data$upper),
                                        colour = "darkorange", inherit.aes = FALSE)
    }
  }
  p <- p + ggplot2::facet_grid(response + .source ~ .covariate, scales = "free",
                              labeller = ggplot2::label_wrap_gen(22))
  p + ggplot2::labs(x = NULL, y = "Response",
                    caption = if (effect != "none") {
                      paste0("Orange: ", effect, " covariate effect (average expected response, ",
                             100 * level, "% credible band); points: fitted values at observed covariates")
                    })
}
