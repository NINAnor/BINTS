# Helpers shared by the plotting functions ------------------------------------

bm_check_model <- function(model) {
  if (!inherits(model, "bm_fit")) stop("`model` must be a bm_fit object.", call. = FALSE)
  invisible(model)
}

bm_check_responses <- function(model, response) {
  bm_check_model(model)
  if (!is.character(response) || !length(response) || anyDuplicated(response) ||
      !all(response %in% model$spec$responses)) {
    stop("`response` must be among the fitted responses: ",
         paste(model$spec$responses, collapse = ", "), ".", call. = FALSE)
  }
  response
}

# Observed values of response r for original-data row indices
bm_obs_for <- function(model, r, rows) {
  model$prep$y[[r]][match(rows, model$prep$kept_rows)]
}

# Column of the modelling data, aligned to original-data row indices
bm_data_col <- function(model, col, rows) {
  model$prep$data[[col]][match(rows, model$prep$kept_rows)]
}

bm_check_facets <- function(n, max_facets, what) {
  if (n > max_facets) {
    stop("`", what, "` has ", n, " levels, more than `max_facets` = ", max_facets,
         ". Subset the data or raise `max_facets`.", call. = FALSE)
  }
}

# Covariate handling --------------------------------------------------------
#
# A single ggplot cannot mix continuous and discrete x scales across facets.
# If all covariates are continuous or all are categorical the natural scale is
# used. For a mix, categorical covariates are drawn at integer positions (in
# level order) and their levels are listed in the facet label.

bm_cov_info <- function(dat, covariates, max_numeric_levels) {
  if (!is.character(covariates) || !length(covariates)) {
    stop("`covariates` must be a character vector of column names.", call. = FALSE)
  }
  miss <- setdiff(covariates, names(dat))
  if (length(miss)) {
    stop("Covariate(s) not found in the modelling data: ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  }
  info <- lapply(stats::setNames(covariates, covariates), function(cv) {
    v <- dat[[cv]]
    u <- unique(v[!is.na(v)])
    cat <- !is.numeric(v) || length(u) <= max_numeric_levels
    lev <- if (!cat) NULL else if (is.factor(v)) levels(droplevels(v)) else {
      if (is.numeric(v)) as.character(sort(u)) else sort(as.character(u))
    }
    list(name = cv, cat = cat, levels = lev)
  })
  cats <- vapply(info, `[[`, logical(1), "cat")
  mode <- if (all(!cats)) "numeric" else if (all(cats)) "factor" else "position"
  list(covariates = info, mode = mode)
}

bm_cov_label <- function(ci, mode) {
  if (mode == "position" && ci$cat) {
    lv <- ci$levels
    if (length(lv) <= 4L) paste0(ci$name, " (", paste(lv, collapse = " | "), ")")
    else paste0(ci$name, " (", length(lv), " levels)")
  } else ci$name
}

bm_x_value <- function(v, ci, mode) {
  if (!ci$cat) return(as.numeric(v))
  if (mode == "factor") return(factor(as.character(v), levels = ci$levels))
  match(as.character(v), ci$levels)
}

# Long data frame with `.covariate` (facet label), `.x`, `.cat` and the
# columns in `keep`, replicated once per covariate.
bm_cov_long <- function(df, cinfo, keep) {
  mode <- cinfo$mode
  parts <- lapply(cinfo$covariates, function(ci) {
    out <- df[keep]
    out$.x <- bm_x_value(df[[ci$name]], ci, mode)
    out$.covariate <- bm_cov_label(ci, mode)
    out$.cat <- ci$cat
    out
  })
  out <- do.call(rbind, parts)
  out$.covariate <- factor(out$.covariate,
                           levels = vapply(cinfo$covariates, bm_cov_label, "", mode = mode))
  if (mode == "factor") {
    all_lev <- unique(unlist(lapply(cinfo$covariates, `[[`, "levels")))
    out$.x <- factor(as.character(out$.x), levels = all_lev)
  }
  rownames(out) <- NULL
  out
}

# Points, boxplots (categorical) and smoother (continuous) for `.y` against `.x`
bm_cov_layers <- function(long, smooth, smooth_se, max_levels, point_alpha = 0.35) {
  layers <- list()
  cont <- long[!long$.cat, , drop = FALSE]
  catg <- long[long$.cat, , drop = FALSE]
  aes_xy <- ggplot2::aes(x = .data$.x, y = .data$.y)
  if (nrow(cont)) {
    layers <- c(layers, list(ggplot2::geom_point(data = cont, mapping = aes_xy,
                                                 alpha = point_alpha, size = 0.9)))
    if (smooth) {
      layers <- c(layers, list(ggplot2::geom_smooth(
        data = cont, mapping = aes_xy, method = "gam", formula = y ~ s(x, bs = "cs"),
        se = smooth_se, colour = "steelblue", fill = "steelblue", alpha = 0.2,
        linewidth = 0.7)))
    }
  }
  if (nrow(catg)) {
    nlev <- tapply(as.character(catg$.x), catg$.covariate, function(v) length(unique(v)))
    many <- catg$.covariate %in% names(nlev)[!is.na(nlev) & nlev > max_levels]
    few <- catg[!many, , drop = FALSE]
    mny <- catg[many, , drop = FALSE]
    if (nrow(few)) {
      layers <- c(layers, list(
        ggplot2::geom_boxplot(data = few, mapping = ggplot2::aes(x = .data$.x, y = .data$.y,
                                                                 group = .data$.x),
                              outlier.shape = NA, fill = NA, colour = "grey30", width = 0.6),
        ggplot2::geom_point(data = few, mapping = aes_xy, alpha = point_alpha, size = 0.9,
                            position = ggplot2::position_jitter(width = 0.15, height = 0,
                                                                seed = 1))))
    }
    if (nrow(mny)) {
      ms <- function(v) {
        m <- mean(v)
        se <- stats::sd(v) / sqrt(length(v))
        data.frame(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
      }
      layers <- c(layers, list(ggplot2::stat_summary(
        data = mny, mapping = aes_xy, fun.data = ms, geom = "pointrange", size = 0.15,
        colour = "grey30")))
    }
  }
  layers
}

# Square axes within a panel so a 1:1 line is meaningful: a blank layer that
# forces x and y to share limits in each panel
bm_square_blank <- function(df, panel_vars, xcol = "observed", ycols = "predicted") {
  key <- do.call(paste, c(df[panel_vars], sep = "\r"))
  parts <- lapply(split(seq_len(nrow(df)), key), function(ii) {
    lim <- range(unlist(df[ii, c(xcol, ycols), drop = FALSE]), na.rm = TRUE)
    out <- df[ii[1], panel_vars, drop = FALSE][c(1, 1), , drop = FALSE]
    out$.lim <- lim
    out
  })
  do.call(rbind, parts)
}
