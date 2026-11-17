# NIMBLE model construction. The model code is assembled from small components
# (intercept, fixed effects, temporal, spatial, random effects, likelihood) as
# language objects, so new components are added rather than model code being
# duplicated for each combination.
#
# Linear predictor for response r, observation i:
#   mu_r[i] = beta0_r + X_r[i, ] %*% beta_r             fixed effects
#           + beta_trend_r * (year - mean year)          linear trend, or
#           + f_cyc_r(season)                            cyclic spline
#           + field_r[site(i)]                           spatial effect
#           + field_spt_r[site(i)] * (year_c | f_cyc)    space-time interaction
#           + sum_e u_{e,r}[id_e[i]]                     random effects
# Truncated normal:  y ~ N(mu, sigma^2) truncated to [0, Inf), via the custom
#                    dtnorm_lower() (see bm_distributions.R)
# Log-normal:        log(y) ~ N(mu, sigma^2)
#
# Sharing (two responses): the reference response (the first) owns the latent
# effect; the second response's version is the deterministic expression
# `lambda * reference effect` (or the same value if `shared_loadings = FALSE`).
# See bm_sharing.R.

bm_response_block <- function(r, prep, spec, priors) {
  S <- function(p) as.name(paste(p, r, sep = "_"))
  nm <- function(p) paste(p, r, sep = "_")
  pr <- priors[[r]]
  fam <- bm_families()[[spec$family[[r]]]]
  sh <- spec$sharing
  if (is.null(sh)) sh <- bm_sharing_inactive()
  ref <- sh$ref
  second <- isTRUE(sh$active) && r == sh$second
  # Loading expression: `lambda * x`, or just `x` for identical sharing
  with_loading <- function(ld, e) if (is.null(ld)) e else bquote(.(ld$sym) * .(e))
  new_loading <- function(stub) {
    if (sh$loadings) bm_loading(stub, r, pr$loading_sd) else NULL
  }

  stmts <- list()
  terms <- list(S("beta0"))
  const <- list()
  monitors <- nm("beta0")

  const[[nm("prior_beta0_mean")]] <- pr$beta0_mean
  const[[nm("prior_beta0_sd")]] <- pr$beta0_sd
  stmts <- c(stmts, list(
    bquote(.(S("beta0")) ~ dnorm(.(S("prior_beta0_mean")),
                                 sd = .(S("prior_beta0_sd"))))
  ))

  add_loading <- function(ld) {
    if (is.null(ld)) return(invisible())
    stmts <<- c(stmts, list(ld$stmt))
    const <<- c(const, ld$const)
    monitors <<- c(monitors, ld$name)
  }

  # Fixed effects
  des <- prep$design[[r]]
  free_cols <- integer()
  if (!is.null(des)) {
    p <- ncol(des$X)
    const[[nm("p")]] <- p
    const[[nm("X")]] <- des$X
    const[[nm("prior_beta_sd")]] <- pr$beta_sd
    shared_cols <- if (second) which(des$scaling$variable %in% sh$fixed) else integer()
    free_cols <- setdiff(seq_len(p), shared_cols)
    if (!length(shared_cols)) {
      stmts <- c(stmts, list(
        bquote(for (k in 1:.(S("p"))) {
          .(S("beta"))[k] ~ dnorm(0, sd = .(S("prior_beta_sd")))
        })
      ))
    } else {
      for (k in free_cols) {
        stmts <- c(stmts, list(
          bquote(.(S("beta"))[.(k)] ~ dnorm(0, sd = .(S("prior_beta_sd"))))
        ))
      }
      for (k in shared_cols) {
        kref <- match(colnames(des$X)[k], colnames(prep$design[[ref]]$X))
        ld <- new_loading(paste0("beta", kref))
        add_loading(ld)
        bref <- as.name(paste("beta", ref, sep = "_"))
        stmts <- c(stmts, list(
          bquote(.(S("beta"))[.(k)] <- .(with_loading(ld, bquote(.(bref)[.(kref)]))))
        ))
      }
    }
    # Explicit sum: inprod() over 1:p fails to compile in NIMBLE when p = 1
    terms <- c(terms, lapply(seq_len(p), function(k) {
      bquote(.(S("beta"))[.(k)] * .(S("X"))[i, .(k)])
    }))
    monitors <- c(monitors, nm("beta"))
  }

  # Structured temporal effect
  trend_derived <- FALSE
  cyc_derived <- FALSE
  if (spec$temporal_effect != "none") {
    tderived <- bm_is_derived(sh, r, "temporal")
    if (spec$temporal_effect == "linear") {
      if (tderived) {
        ld <- new_loading("trend")
        add_loading(ld)
        bt_ref <- as.name(paste("beta_trend", ref, sep = "_"))
        stmts <- c(stmts, list(
          bquote(.(S("beta_trend")) <- .(with_loading(ld, bt_ref)))
        ))
        trend_derived <- TRUE
      } else {
        const[[nm("prior_trend_sd")]] <- pr$trend_sd
        stmts <- c(stmts, list(
          bquote(.(S("beta_trend")) ~ dnorm(0, sd = .(S("prior_trend_sd"))))
        ))
      }
      terms <- c(terms, list(bquote(.(S("beta_trend")) * year_c[i])))
      monitors <- c(monitors, nm("beta_trend"))
    } else {
      if (tderived) {
        ld <- new_loading("cyc")
        add_loading(ld)
        rw <- as.name(paste("w_cyc", ref, sep = "_"))
        rs <- as.name(paste("sd_cyc", ref, sep = "_"))
        stmts <- c(stmts, list(
          bquote(.(S("sd_cyc")) <- .(if (is.null(ld)) rs else bquote(abs(.(ld$sym)) * .(rs)))),
          bquote(.(S("w_cyc"))[1:K_cyc] <- .(with_loading(ld, bquote(.(rw)[1:K_cyc]))))
        ))
        cyc_derived <- TRUE
      } else {
        const[[nm("prior_cyc_tau")]] <- 1 / pr$smooth_scale^2
        stmts <- c(stmts, list(
          bquote(.(S("sd_cyc")) ~ T(dt(0, tau = .(S("prior_cyc_tau")), df = 3), 0, Inf)),
          bquote(for (k in 1:K_cyc) .(S("w_cyc"))[k] ~ dnorm(0, sd = .(S("sd_cyc"))))
        ))
      }
      stmts <- c(stmts, list(
        bquote(.(S("f_cyc"))[1:n] <- Bcyc[1:n, 1:K_cyc] %*% .(S("w_cyc"))[1:K_cyc])
      ))
      terms <- c(terms, list(bquote(.(S("f_cyc"))[i])))
      monitors <- c(monitors, nm("sd_cyc"), nm("w_cyc"))
    }
  }

  # Structured spatial effect (field defined at sites; see bm_spatial.R)
  fields <- character()
  nngp_fields <- list()
  share_to <- if (isTRUE(sh$active) && r == ref && isTRUE(sh$spatial)) sh$second
  sp_derived <- bm_is_derived(sh, r, "spatial")
  derived_arg <- if (sp_derived) {
    list(ref = ref, loading_sd = if (sh$loadings) pr$loading_sd)
  }
  if (spec$spatial_effect != "none") {
    sp_blk <- bm_spatial_block(r, prep, pr$spatial_scale, pr$beta_sd, "sp",
                               derived = derived_arg, share_to = share_to)
    nngp_fields <- c(nngp_fields, sp_blk$nngp)
    stmts <- c(stmts, sp_blk$stmts)
    terms <- c(terms, sp_blk$terms)
    const <- c(const, sp_blk$constants)
    monitors <- c(monitors, sp_blk$monitors)
    fields <- "sp"
  }

  # Space-time interaction: the spatial pattern of the temporal term varies
  # over space (a spatially varying slope or seasonal amplitude)
  if (isTRUE(spec$spatiotemporal_interaction)) {
    mult <- if (spec$temporal_effect == "linear") quote(year_c[i]) else {
      bquote(.(S("f_cyc"))[i])
    }
    st_blk <- bm_spatial_block(r, prep, pr$interaction_scale, pr$interaction_scale,
                               "spt", mult, derived = derived_arg, share_to = share_to)
    nngp_fields <- c(nngp_fields, st_blk$nngp)
    stmts <- c(stmts, st_blk$stmts)
    terms <- c(terms, st_blk$terms)
    const <- c(const, st_blk$constants)
    monitors <- c(monitors, st_blk$monitors)
    fields <- c(fields, "spt")
  }

  # Random effects (centred; their means are absorbed by the intercept)
  if (length(spec$random_effects)) {
    const[[nm("prior_re_tau")]] <- 1 / pr$sd_re_scale^2
  }
  year_derived <- FALSE
  for (e in spec$random_effects) {
    u <- S(paste0("u_", e))
    s <- S(paste0("sd_", e))
    N <- as.name(paste0("n_", e))
    id <- as.name(paste0(e, "_id"))
    u_ref <- as.name(paste("u", e, ref, sep = "_"))
    s_ref <- as.name(paste("sd", e, ref, sep = "_"))
    if (bm_is_derived(sh, r, "random", e)) {
      ld <- new_loading(paste0("re_", e))
      add_loading(ld)
      stmts <- c(stmts, list(
        bquote(for (j in 1:.(N)) .(u)[j] <- .(with_loading(ld, bquote(.(u_ref)[j])))),
        bquote(.(s) <- .(if (is.null(ld)) s_ref else bquote(abs(.(ld$sym)) * .(s_ref))))
      ))
      if (e == "year") year_derived <- TRUE
    } else if (second && e %in% sh$correlated) {
      rho <- as.name(paste("rho", e, r, sep = "_"))
      stmts <- c(stmts, list(
        bquote(for (j in 1:.(N)) {
          .(u)[j] ~ dnorm(.(rho) * .(s) / .(s_ref) * .(u_ref)[j],
                          sd = .(s) * sqrt(1 - .(rho) * .(rho)))
        }),
        bquote(.(s) ~ T(dt(0, tau = .(S("prior_re_tau")), df = 3), 0, Inf)),
        bquote(.(rho) ~ dunif(-1, 1))
      ))
      monitors <- c(monitors, as.character(rho))
    } else {
      stmts <- c(stmts, list(
        bquote(for (j in 1:.(N)) .(u)[j] ~ dnorm(0, sd = .(s))),
        bquote(.(s) ~ T(dt(0, tau = .(S("prior_re_tau")), df = 3), 0, Inf))
      ))
    }
    terms <- c(terms, list(bquote(.(u)[.(id)[i]])))
    monitors <- c(monitors, nm(paste0("sd_", e)), nm(paste0("u_", e)))
  }

  # Likelihood. The family's second parameter is the residual sd `sigma` or, for
  # the negative binomial, the size `phi`.
  p2 <- fam$par2
  par2_stmt <- NULL
  if (is.null(p2)) {
    # no second parameter (Poisson)
  } else if (identical(fam$par2_prior, "gamma")) {
    const[[nm("prior_disp_shape")]] <- pr$disp_shape
    const[[nm("prior_disp_rate")]] <- pr$disp_rate
    par2_stmt <- bquote(.(S(p2)) ~ dgamma(.(S("prior_disp_shape")), .(S("prior_disp_rate"))))
  } else {
    const[[nm("prior_sigma_tau")]] <- 1 / pr$sigma_scale^2
    par2_stmt <- bquote(.(S(p2)) ~ T(dt(0, tau = .(S("prior_sigma_tau")), df = 3), 0, Inf))
  }
  mu_expr <- Reduce(function(a, b) call("+", a, b), terms)
  lik <- fam$likelihood(S("y"), S("mu"), if (is.null(p2)) NULL else S(p2))
  lik_stmts <- if (is.call(lik) && identical(lik[[1]], as.name("{"))) as.list(lik)[-1] else list(lik)
  loop_body <- as.call(c(as.name("{"), c(list(bquote(.(S("mu"))[i] <- .(mu_expr))), lik_stmts)))
  stmts <- c(stmts, if (!is.null(par2_stmt)) list(par2_stmt),
             list(call("for", quote(i), quote(1:n), loop_body)))
  if (!is.null(p2)) monitors <- c(monitors, nm(p2))

  # --- Samplers -------------------------------------------------------------
  # Block samplers for strongly correlated low-dimensional parameters: the
  # intercept, fixed/trend coefficients and the (few) year effects are
  # confounded through their sum; cyclic/spline coefficients are mutually
  # correlated. Only parameters this response owns (non-derived) are used.
  # With an NNGP the intercept (and, with an interaction, the trend) is moved by
  # shift samplers instead; their posteriors are wide, which would distort the
  # block proposal. When the reference response's latent is shared, the shift
  # also moves the second response's (derived) effect, so the second response's
  # free node is shifted to compensate (see bm_distributions.R).
  second_name <- if (isTRUE(sh$active)) sh$second
  ref_owns <- isTRUE(sh$active) && r == ref
  lam_of <- function(stub) {
    if (sh$loadings) paste("lambda", stub, second_name, sep = "_")
  }
  nngp <- spec$spatial_effect == "autocorrelation"
  inter <- isTRUE(spec$spatiotemporal_interaction)
  trend_shift <- nngp && inter && spec$temporal_effect == "linear" &&
    !sp_derived && !isTRUE(sh$temporal)
  shifts <- list()
  # the whole field is one vector node, addressed by its variable name
  site_nodes <- function(tg) nm(paste0("w_", tg))
  mk_shift <- function(target, nodes, scale, shared_flag, stub, target2) {
    s <- list(target = target, nodes = nodes, scale = scale)
    if (ref_owns && shared_flag) {
      s$target2 <- target2
      s$lamNode <- lam_of(stub)
    }
    s
  }
  if (nngp && !sp_derived) {
    shifts <- c(shifts, list(mk_shift(nm("beta0"), site_nodes("sp"),
                                      0.2 * pr$sigma_scale, sh$spatial, "sp",
                                      paste("beta0", second_name, sep = "_"))))
    if ("year" %in% spec$random_effects && !year_derived) {
      shifts <- c(shifts, list(mk_shift(
        nm("beta0"),
        paste0(nm("u_year"), "[", seq_along(prep$index$year$levels), "]"),
        0.2 * pr$sigma_scale, "year" %in% sh$random, "re_year",
        paste("beta0", second_name, sep = "_"))))
    }
  }
  if (trend_shift) {
    shifts <- c(shifts, list(mk_shift(nm("beta_trend"), site_nodes("spt"),
                                      0.05 * pr$trend_sd, sh$spatial, "spt",
                                      paste("beta_trend", second_name, sep = "_"))))
  }

  blocks <- list()
  core <- if (nngp) character() else nm("beta0")
  # (paste0() with an empty index would give the whole, partly derived vector)
  if (!is.null(des) && length(free_cols)) {
    core <- c(core, paste0(nm("beta"), "[", free_cols, "]"))
  }
  if (spec$temporal_effect == "linear" && !trend_derived && !trend_shift) {
    core <- c(core, nm("beta_trend"))
  }
  if ("year" %in% spec$random_effects && !year_derived &&
      !(second && "year" %in% sh$correlated)) {
    core <- c(core, paste0(nm("u_year"), "[", seq_along(prep$index$year$levels), "]"))
  }
  if (length(core) > 1L) blocks <- c(blocks, list(core))
  if (spec$spatial_effect == "coordinates" && !sp_derived) {
    for (tg in fields) {
      blocks <- c(blocks, list(paste0(nm(paste0("v_", tg)), "[", seq_len(prep$spatial$tp$K), "]")))
    }
  }
  if (spec$temporal_effect == "cyclic" && !cyc_derived) {
    blocks <- c(blocks, list(paste0(nm("w_cyc"), "[", seq_len(prep$temporal$cyclic$K), "]")))
  }

  list(stmts = stmts, constants = const, data = stats::setNames(list(prep$y[[r]]), nm("y")),
       monitors = monitors, blocks = blocks, shifts = shifts, nngp = nngp_fields)
}

# Parameters (with dimensions) that need initial values, per response. Derived
# (deterministic) nodes get none; loadings start at half the scale ratio.
bm_make_inits <- function(prep, spec, priors, seed) {
  sh <- spec$sharing
  if (is.null(sh)) sh <- bm_sharing_inactive()
  ref <- sh$ref
  bm_with_seed(seed, {
    inits <- list()
    for (r in spec$responses) {
      pr <- priors[[r]]
      sc <- pr$sigma_scale
      second <- isTRUE(sh$active) && r == sh$second
      lam_init <- 0.5 * pr$loading_sd / 5
      inits[[paste0("beta0_", r)]] <- pr$beta0_mean + stats::rnorm(1, 0, 0.1 * sc)
      des <- prep$design[[r]]
      if (!is.null(des)) {
        shared_cols <- if (second) which(des$scaling$variable %in% sh$fixed) else integer()
        free <- setdiff(seq_len(ncol(des$X)), shared_cols)
        b <- rep(NA_real_, ncol(des$X))
        b[free] <- stats::rnorm(length(free), 0, 0.05 * sc)
        inits[[paste0("beta_", r)]] <- b
        if (sh$loadings) {
          for (k in shared_cols) {
            kref <- match(colnames(des$X)[k], colnames(prep$design[[ref]]$X))
            inits[[paste0("lambda_beta", kref, "_", r)]] <- lam_init
          }
        }
      }
      tder <- bm_is_derived(sh, r, "temporal")
      if (spec$temporal_effect == "linear") {
        if (tder) {
          if (sh$loadings) inits[[paste0("lambda_trend_", r)]] <- lam_init
        } else inits[[paste0("beta_trend_", r)]] <- 0
      } else if (spec$temporal_effect == "cyclic") {
        if (tder) {
          if (sh$loadings) inits[[paste0("lambda_cyc_", r)]] <- lam_init
        } else {
          inits[[paste0("sd_cyc_", r)]] <- sc * stats::runif(1, 0.1, 0.4)
          inits[[paste0("w_cyc_", r)]] <-
            stats::rnorm(prep$temporal$cyclic$K, 0, 0.1 * sc)
        }
      }
      tags <- c(if (spec$spatial_effect != "none") "sp",
                if (isTRUE(spec$spatiotemporal_interaction)) "spt")
      for (tg in tags) {
        if (bm_is_derived(sh, r, "spatial")) {
          if (sh$loadings) inits[[paste0("lambda_", tg, "_", r)]] <- lam_init
          next
        }
        # scale of the field: the response for the main field, the interaction
        # scale (slope per year or amplitude multiplier) for the interaction
        fs <- if (tg == "sp") sc else priors[[r]]$interaction_scale
        inits[[paste0("sd_", tg, "_", r)]] <- fs * stats::runif(1, 0.1, 0.4)
        if (spec$spatial_effect == "autocorrelation") {
          ng <- prep$spatial$nngp
          inits[[paste0("phi_idx_", tg, "_", r)]] <-
            sample(ceiling(ng$G / 3):ceiling(2 * ng$G / 3), 1L)
          inits[[paste0("w_", tg, "_", r)]] <- stats::rnorm(nrow(ng$xy), 0, 0.05 * fs)
        } else {
          inits[[paste0("v_", tg, "_", r)]] <- stats::rnorm(prep$spatial$tp$K, 0, 0.05 * fs)
        }
      }
      for (e in spec$random_effects) {
        if (bm_is_derived(sh, r, "random", e)) {
          if (sh$loadings) inits[[paste0("lambda_re_", e, "_", r)]] <- lam_init
          next
        }
        inits[[paste0("sd_", e, "_", r)]] <- sc * stats::runif(1, 0.2, 0.6)
        inits[[paste0("u_", e, "_", r)]] <-
          stats::rnorm(length(prep$index[[e]]$levels), 0, 0.05 * sc)
        if (second && e %in% sh$correlated) inits[[paste0("rho_", e, "_", r)]] <- 0
      }
      fam_r <- bm_families()[[spec$family[[r]]]]
      if (!is.null(fam_r$par2)) {
        inits[[paste0(fam_r$par2, "_", r)]] <- fam_r$par2_init(sc, prep$y[[r]])
      }
    }
    inits
  })
}

bm_build_model <- function(prep, spec, priors) {
  if (is.null(spec$temporal_effect)) spec$temporal_effect <- "none"
  if (is.null(spec$spatial_effect)) spec$spatial_effect <- "none"
  blocks <- lapply(spec$responses, bm_response_block, prep = prep, spec = spec,
                   priors = priors)
  stmts <- unlist(lapply(blocks, `[[`, "stmts"), recursive = FALSE)
  body <- as.call(c(as.name("{"), stmts))
  code <- eval(substitute(nimble::nimbleCode(CODE), list(CODE = body)))

  constants <- c(
    list(n = prep$n),
    stats::setNames(lapply(prep$index, function(x) length(x$levels)),
                    paste0("n_", names(prep$index))),
    stats::setNames(lapply(prep$index, `[[`, "id"),
                    paste0(names(prep$index), "_id")),
    if (spec$temporal_effect == "linear") list(year_c = prep$temporal$year_c),
    if (spec$spatial_effect != "none") bm_spatial_constants(prep),
    if (spec$temporal_effect == "cyclic") {
      list(Bcyc = prep$temporal$cyclic$basis, K_cyc = prep$temporal$cyclic$K)
    },
    unlist(lapply(blocks, `[[`, "constants"), recursive = FALSE)
  )
  list(
    code = code,
    constants = constants,
    data = c(unlist(lapply(blocks, `[[`, "data"), recursive = FALSE),
             if (spec$spatial_effect == "autocorrelation") {
               list(Bflat = prep$spatial$nngp$Bflat, sqrtF = prep$spatial$nngp$sqrtF)
             }),
    monitors = unlist(lapply(blocks, `[[`, "monitors")),
    blocks = unlist(lapply(blocks, `[[`, "blocks"), recursive = FALSE),
    shifts = unlist(lapply(blocks, `[[`, "shifts"), recursive = FALSE),
    nngp = if (spec$spatial_effect == "autocorrelation") {
      ng <- prep$spatial$nngp
      list(fields = unlist(lapply(blocks, `[[`, "nngp"), recursive = FALSE),
           L = nrow(ng$xy), site_id = prep$spatial$site_id, nbr = ng$nbr,
           Bflat = ng$Bflat, sqrtF = ng$sqrtF, child = ng$child, nChild = ng$nChild)
    }
  )
}
