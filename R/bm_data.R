#' Default column mapping for biomass-richness models
#'
#' Returns the mapping between the model's roles and the column names in the
#' data. Pass a modified copy to the `columns` argument of [fit_model()] to use
#' other column names; nothing in the package refers to column names directly.
#'
#' The `replication` column is converted to a replicate index using its
#' trailing number (the text after the final non-digit, e.g. `"..._sampling-03"`
#' gives 3), so replicate occasions are crossed with locality and year. If the
#' column has no trailing number, its unique values are used as levels.
#'
#' @return A named list of column names.
#' @export
#' @examples
#' bm_default_columns()
#' @md
bm_default_columns <- function() {
  list(
    year = "year",
    locality = "locality",
    replication = "sampling_name",
    season = "start_date_julian",
    longitude = "longitude",
    latitude = "latitude",
    biomass = "sum_wet_weight",
    richness = "no_species"
  )
}

bm_make_index <- function(x, effect) {
  if (effect == "replication") {
    num <- suppressWarnings(
      as.integer(sub(".*[^0-9]([0-9]+)$", "\\1", as.character(x)))
    )
    if (!anyNA(num)) x <- num
  }
  f <- factor(x)
  list(id = as.integer(f), levels = levels(f))
}

# Build the design matrix for one response. Numeric columns with more than two
# distinct values are centred and scaled; dummy columns are left as 0/1.
bm_design_matrix <- function(dat, vars) {
  if (!length(vars)) return(NULL)
  f <- stats::as.formula(paste("~", paste0("`", vars, "`", collapse = " + ")))
  mm_full <- stats::model.matrix(f, data = dat)
  assign_idx <- attr(mm_full, "assign")[-1]
  mm <- mm_full[, -1, drop = FALSE]
  vars_by_col <- attr(stats::terms(f), "term.labels")[assign_idx]
  vars_by_col <- gsub("`", "", vars_by_col, fixed = TRUE)
  scaling <- data.frame(
    term = colnames(mm), variable = vars_by_col, center = 0, scale = 1,
    row.names = NULL
  )
  for (j in seq_len(ncol(mm))) {
    if (length(unique(mm[, j])) > 2L) {
      scaling$center[j] <- mean(mm[, j])
      scaling$scale[j] <- stats::sd(mm[, j])
      if (scaling$scale[j] == 0) {
        stop("Fixed effect '", colnames(mm)[j], "' is constant.", call. = FALSE)
      }
      mm[, j] <- (mm[, j] - scaling$center[j]) / scaling$scale[j]
    }
  }
  list(X = mm, scaling = scaling)
}

# Validate columns, drop incomplete rows, and construct indices and design
# matrices used by the model builder.
bm_prepare_data <- function(data, responses, fixed_effects, random_effects,
                            family, columns, temporal_effect = "none",
                            temporal_options = NULL, spatial_effect = "none",
                            spatial_options = NULL) {
  fam <- bm_families()
  resp_cols <- unlist(columns[responses])
  re_cols <- unlist(columns[random_effects])
  fx_cols <- unique(unlist(fixed_effects))
  tm_cols <- switch(temporal_effect, none = NULL, linear = columns$year,
                    cyclic = columns$season)
  if (temporal_effect != "none" && is.null(tm_cols)) {
    stop("`columns` has no entry for the ", temporal_effect,
         " temporal variable ('", if (temporal_effect == "linear") "year" else "season",
         "').", call. = FALSE)
  }
  sp_cols <- if (spatial_effect != "none") {
    c(columns$locality, columns$longitude, columns$latitude)
  }
  needed <- unique(c(resp_cols, re_cols, fx_cols, tm_cols, sp_cols))

  miss <- setdiff(needed, names(data))
  if (length(miss)) {
    stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  for (r in responses) {
    if (!is.numeric(data[[columns[[r]]]])) {
      stop("Response column '", columns[[r]], "' must be numeric.", call. = FALSE)
    }
  }

  keep <- stats::complete.cases(as.data.frame(data)[, needed, drop = FALSE])
  if (!all(keep)) {
    message(sum(!keep), " row(s) with missing values in model variables were ",
            "dropped.")
  }
  if (sum(keep) < 10L) stop("Fewer than 10 complete rows available.", call. = FALSE)
  dat <- as.data.frame(data)[keep, , drop = FALSE]
  rownames(dat) <- NULL

  y <- list()
  for (r in responses) {
    yr <- dat[[columns[[r]]]]
    f <- fam[[family[[r]]]]
    if (!f$valid(yr)) {
      stop("Response '", r, "' (column '", columns[[r]], "') must contain ",
           f$valid_msg, " for the ", f$label, " family.",
           if (family[[r]] == "lognormal" && any(yr == 0)) {
             paste0(" It contains ", sum(yr == 0), " zero value(s); use ",
                    "family = 'truncnorm' or remove/offset the zeros.")
           } else "", call. = FALSE)
    }
    y[[r]] <- yr
  }

  index <- lapply(stats::setNames(random_effects, random_effects), function(e) {
    bm_make_index(dat[[columns[[e]]]], e)
  })

  design <- lapply(fixed_effects[responses], function(v) bm_design_matrix(dat, v))

  temporal_options <- bm_resolve_temporal_options(temporal_options)
  temporal <- bm_prepare_temporal(dat, columns, temporal_effect, temporal_options)

  spatial <- bm_prepare_spatial(dat, columns, spatial_effect, spatial_options)

  list(
    n = nrow(dat), y = y, index = index, design = design, data = dat,
    temporal = temporal, spatial = spatial, kept_rows = which(keep), n_dropped = sum(!keep)
  )
}
