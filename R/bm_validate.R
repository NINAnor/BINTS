# Argument validation for fit_model(). Features scheduled for later
# implementation stages are accepted in the interface but rejected here with an
# informative message.

bm_validate_args <- function(data, responses, fixed_effects, random_effects,
                             temporal_effect, spatial_effect,
                             spatiotemporal_interaction, shared, coordinates,
                             columns, n_chains, n_iter, n_burnin, n_thin,
                             seed, chain_seeds, parallel) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (inherits(data, "sf")) {
    stop("`data` must be a plain data frame; supply `coordinates` as column names.",
         call. = FALSE)
  }

  ok_resp <- c("biomass", "richness")
  if (!is.character(responses) || !length(responses) ||
      !all(responses %in% ok_resp) || anyDuplicated(responses)) {
    stop("`responses` must be one or both of: ",
         paste(ok_resp, collapse = ", "), ".", call. = FALSE)
  }

  ok_re <- c("year", "locality", "replication")
  if (!is.null(random_effects) &&
      (!is.character(random_effects) || !all(random_effects %in% ok_re))) {
    stop("`random_effects` must be NULL or a subset of: ",
         paste(ok_re, collapse = ", "), ".", call. = FALSE)
  }

  if (!is.list(columns) || is.null(names(columns))) {
    stop("`columns` must be a named list; see bm_default_columns().",
         call. = FALSE)
  }
  if (spatial_effect != "none") {
    miss_xy <- setdiff(unlist(columns[c("locality", "longitude", "latitude")]),
                       names(data))
    if (length(miss_xy)) {
      stop("Spatial effects need column(s) not found in `data`: ",
           paste(miss_xy, collapse = ", "), ". Use `coordinates` or `columns`.",
           call. = FALSE)
    }
  }
  needed_keys <- c(responses, random_effects)
  no_map <- setdiff(needed_keys, names(columns))
  if (length(no_map)) {
    stop("`columns` has no entry for: ", paste(no_map, collapse = ", "), ".",
         call. = FALSE)
  }

  if (!is.character(temporal_effect) || length(temporal_effect) != 1L ||
      !temporal_effect %in% c("none", "linear", "cyclic")) {
    stop("`temporal_effect` must be 'none', 'linear' or 'cyclic'.", call. = FALSE)
  }
  if (!is.character(spatial_effect) || length(spatial_effect) != 1L ||
      !spatial_effect %in% c("none", "autocorrelation", "coordinates")) {
    stop("`spatial_effect` must be 'none', 'autocorrelation' or 'coordinates'.",
         call. = FALSE)
  }
  if (!is.logical(spatiotemporal_interaction) ||
      length(spatiotemporal_interaction) != 1L || is.na(spatiotemporal_interaction)) {
    stop("`spatiotemporal_interaction` must be TRUE or FALSE.", call. = FALSE)
  }

  not_yet <- character()
  if (spatiotemporal_interaction && (temporal_effect == "none" || spatial_effect == "none")) {
    stop("`spatiotemporal_interaction = TRUE` needs both a temporal effect ",
         "('linear' or 'cyclic') and a spatial effect.", call. = FALSE)
  }
  if (length(not_yet)) {
    stop("Not yet implemented: ", paste(not_yet, collapse = ", "),
         " (planned in later development stages).", call. = FALSE)
  }
  is_set <- vapply(shared, function(s) !isFALSE(s), logical(1))
  if (any(is_set)) {
    if (length(responses) < 2L) {
      warning("Sharing arguments (", paste(names(shared)[is_set], collapse = ", "),
              ") are ignored for a single response.", call. = FALSE)
    }
  }

  one_int <- function(x, nm, min = 1) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min ||
        x != round(x)) {
      stop("`", nm, "` must be a single integer >= ", min, ".", call. = FALSE)
    }
  }
  one_int(n_chains, "n_chains")
  one_int(n_iter, "n_iter")
  one_int(n_burnin, "n_burnin", 0)
  one_int(n_thin, "n_thin")
  if (n_burnin >= n_iter) stop("`n_burnin` must be smaller than `n_iter`.", call. = FALSE)
  if ((n_iter - n_burnin) %/% n_thin < 10) {
    stop("Fewer than 10 saved iterations per chain; increase `n_iter` or ",
         "reduce `n_burnin`/`n_thin`.", call. = FALSE)
  }
  if (!is.null(seed)) one_int(seed, "seed", 0)
  if (!is.null(chain_seeds)) {
    if (!is.numeric(chain_seeds) || length(chain_seeds) != n_chains ||
        anyNA(chain_seeds) || anyDuplicated(chain_seeds)) {
      stop("`chain_seeds` must be ", n_chains, " unique integers (one per chain).",
           call. = FALSE)
    }
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

# Normalise `fixed_effects` into a named list with one character vector per
# response. A plain character vector applies to every response.
bm_resolve_fixed <- function(fixed_effects, responses) {
  if (is.null(fixed_effects)) {
    return(stats::setNames(rep(list(character()), length(responses)), responses))
  }
  if (is.character(fixed_effects)) {
    return(stats::setNames(rep(list(fixed_effects), length(responses)), responses))
  }
  if (is.list(fixed_effects) && !is.null(names(fixed_effects)) &&
      all(names(fixed_effects) %in% responses)) {
    out <- stats::setNames(rep(list(character()), length(responses)), responses)
    out[names(fixed_effects)] <- lapply(fixed_effects, as.character)
    return(out)
  }
  stop("`fixed_effects` must be NULL, a character vector, or a list named by ",
       "response (e.g. list(biomass = 'mean_temp', richness = NULL)).",
       call. = FALSE)
}
