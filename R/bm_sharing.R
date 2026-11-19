# Sharing of effects between two responses.
#
# Three distinct ideas are supported (see ?fit_model):
#  1. Shared parameter (`shared_loadings = FALSE`): both responses use the very
#     same coefficient / random-effect values. Only sensible if the responses
#     are on comparable scales.
#  2. Shared latent process (`shared_loadings = TRUE`, the default): one latent
#     effect `z` whose effect on response r is `lambda_r * z`. The first
#     response in `responses` is the reference and owns `z` (lambda = 1, which
#     fixes sign and scale); the second response gets a free loading, so the
#     two can respond with different magnitude, or opposite sign, to the same
#     process. In the model code the second response's effect is a
#     *deterministic* function (`lambda * reference effect`).
#  3. Correlated response-specific effects (`correlated_random_effects`): each
#     response has its own random effect, with a correlation `rho` between
#     them (a bivariate normal).

bm_sharing_inactive <- function() {
  list(active = FALSE, ref = NULL, second = NULL, loadings = TRUE,
       fixed = character(), random = character(), spatial = FALSE,
       temporal = FALSE, correlated = character())
}

bm_flag_or_names <- function(x, nm, universe, what) {
  if (is.null(x) || isFALSE(x)) return(character())
  if (isTRUE(x)) {
    if (!length(universe)) {
      stop("`", nm, " = TRUE` but there is no ", what, " to share.", call. = FALSE)
    }
    return(universe)
  }
  if (!is.character(x)) {
    stop("`", nm, "` must be FALSE, TRUE or a character vector of ", what,
         " names.", call. = FALSE)
  }
  bad <- setdiff(x, universe)
  if (length(bad)) {
    stop("`", nm, "`: ", paste(bad, collapse = ", "), " not among the available ",
         what, "s (", if (length(universe)) paste(universe, collapse = ", ") else "none",
         ").", call. = FALSE)
  }
  x
}

bm_flag <- function(x, nm, available, what) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", nm, "` must be TRUE or FALSE.", call. = FALSE)
  }
  if (x && !available) {
    stop("`", nm, " = TRUE` needs a ", what, ".", call. = FALSE)
  }
  x
}

bm_resolve_sharing <- function(responses, fixed, random_effects, temporal_effect,
                               spatial_effect, shared_fixed_effects,
                               shared_random_effects, shared_spatial_effect,
                               shared_temporal_effect, correlated_random_effects,
                               shared_loadings) {
  if (length(responses) < 2L) return(bm_sharing_inactive())
  if (!is.logical(shared_loadings) || length(shared_loadings) != 1L ||
      is.na(shared_loadings)) {
    stop("`shared_loadings` must be TRUE or FALSE.", call. = FALSE)
  }
  common <- intersect(fixed[[responses[1]]], fixed[[responses[2]]])
  sf <- bm_flag_or_names(shared_fixed_effects, "shared_fixed_effects", common,
                         "fixed effect common to both responses")
  sr <- bm_flag_or_names(shared_random_effects, "shared_random_effects",
                         random_effects, "random effect")
  cr <- bm_flag_or_names(correlated_random_effects, "correlated_random_effects",
                         random_effects, "random effect")
  both <- intersect(sr, cr)
  if (length(both)) {
    stop("Random effect(s) ", paste(both, collapse = ", "), " cannot be both ",
         "shared and correlated.", call. = FALSE)
  }
  list(
    active = TRUE, ref = responses[1], second = responses[2],
    loadings = shared_loadings, fixed = sf, random = sr,
    spatial = bm_flag(shared_spatial_effect, "shared_spatial_effect",
                      spatial_effect != "none", "spatial effect"),
    temporal = bm_flag(shared_temporal_effect, "shared_temporal_effect",
                       temporal_effect != "none", "temporal effect"),
    correlated = cr
  )
}

bm_any_shared <- function(sh) {
  sh$active && (length(sh$fixed) || length(sh$random) || sh$spatial ||
                  sh$temporal || length(sh$correlated))
}

# Is the effect of response `r` derived from the reference response?
bm_is_derived <- function(sh, r, kind, name = NULL) {
  if (is.null(sh) || !isTRUE(sh$active) || r != sh$second) return(FALSE)
  switch(kind,
         fixed = name %in% sh$fixed,
         random = name %in% sh$random,
         spatial = sh$spatial,
         temporal = sh$temporal,
         FALSE)
}

# Loading node `lambda_<stub>_<r> ~ Normal(0, sd)`
bm_loading <- function(stub, r, sd) {
  node <- paste("lambda", stub, r, sep = "_")
  cn <- paste("prior_lsd", stub, r, sep = "_")
  list(sym = as.name(node), name = node,
       stmt = bquote(.(as.name(node)) ~ dnorm(0, sd = .(as.name(cn)))),
       const = stats::setNames(list(sd), cn))
}
