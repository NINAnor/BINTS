# MCMC configuration, execution and seeds ------------------------------------

# Evaluate `expr` with a temporary RNG seed, restoring the caller's RNG state.
bm_with_seed <- function(seed, expr) {
  has_old <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_old) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (has_old) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

# Master seed (recorded) -> one unique seed per chain, unless supplied
bm_chain_seeds <- function(seed, chain_seeds, n_chains) {
  if (!is.null(chain_seeds)) {
    return(list(master = NA_integer_, chains = as.integer(chain_seeds)))
  }
  if (is.null(seed)) seed <- sample.int(1e6, 1L)
  seeds <- bm_with_seed(seed, sample.int(.Machine$integer.max, n_chains))
  list(master = as.integer(seed), chains = seeds)
}

# Build, compile and run MCMC for one or more chains in a single process.
# Defined to depend only on nimble:: so it can run on PSOCK workers without the
# package installed (its environment is reset to globalenv() before export).
bm_run_chains <- function(inits_list, seeds, code, constants, data, monitors,
                          n_iter, n_burnin, n_thin, dist_code, blocks = NULL,
                          shifts = NULL, nngp = NULL) {
  # NIMBLE's internals look up some of its own functions on the search path
  if (!"package:nimble" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("nimble"))
    on.exit(detach("package:nimble", character.only = TRUE), add = TRUE)
  }
  # NIMBLE only finds user-defined distributions in the global environment;
  # define them temporarily and remove them on exit.
  eval(dist_code, globalenv())
  on.exit(rm(list = c("dtnorm_lower", "rtnorm_lower", "sampler_shift", "dnngp", "rnngp",
                      "siteCalc_virtual", "siteCalc", "sampler_nngp"),
             envir = globalenv()), add = TRUE)
  model <- nimble::nimbleModel(code, constants = constants, data = data,
                               inits = inits_list[[1]], check = FALSE)
  cmodel <- nimble::compileNimble(model)
  conf <- nimble::configureMCMC(model, monitors = monitors, print = FALSE)
  stoch <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  for (b in blocks) {
    # only free (stochastic, unobserved) nodes may be sampled in a block
    b <- b[b %in% stoch]
    if (length(b) < 2L) next
    conf$removeSamplers(b)
    conf$addSampler(target = b, type = "RW_block")
  }
  # NNGP fields: replace the default sampler of each vector node
  for (f in nngp$fields) {
    rows <- split(seq_along(nngp$site_id), factor(nngp$site_id, levels = seq_len(nngp$L)))
    site_nodes <- lapply(seq_len(nngp$L), function(l) {
      start <- paste0(f$mu, "[", rows[[l]], "]")
      if (!is.null(f$derived)) start <- c(paste0(f$derived, "[", l, "]"), start)
      model$getDependencies(start, self = TRUE)
    })
    tgt <- model$expandNodeNames(f$target)
    conf$removeSamplers(tgt)
    conf$addSampler(
      target = tgt, type = sampler_nngp,
      control = list(L = nngp$L, nbr = nngp$nbr, Bflat = nngp$Bflat, sqrtF = nngp$sqrtF,
                     child = nngp$child, nChild = nngp$nChild, sdNode = f$sd,
                     idxNode = f$idx, siteNodes = site_nodes, scale0 = f$scale0)
    )
  }
  for (s in shifts) {
    t2 <- if (is.null(s$target2)) s$target else s$target2
    lam <- if (is.null(s$lamNode)) s$target else s$lamNode
    conf$addSampler(
      target = s$target, type = sampler_shift,
      control = list(shiftNodes = s$nodes, scale = s$scale, target2 = t2,
                     lamNode = lam, hasT2 = !is.null(s$target2),
                     hasLam = !is.null(s$lamNode))
    )
  }
  mcmc <- nimble::buildMCMC(conf)
  cmcmc <- nimble::compileNimble(mcmc, project = model)
  lapply(seq_along(seeds), function(k) {
    nimble::runMCMC(cmcmc, niter = n_iter, nburnin = n_burnin, thin = n_thin,
                    inits = inits_list[[k]], setSeed = seeds[k],
                    progressBar = FALSE, samplesAsCodaMCMC = FALSE)
  })
}

bm_run_mcmc <- function(built, inits, seeds, n_iter, n_burnin, n_thin, parallel) {
  n_chains <- length(seeds)
  worker <- bm_run_chains
  environment(worker) <- globalenv()
  args <- list(code = built$code, constants = built$constants, data = built$data,
               monitors = built$monitors, n_iter = n_iter, n_burnin = n_burnin,
               n_thin = n_thin, dist_code = bm_dist_code(),
               blocks = built$blocks, shifts = built$shifts, nngp = built$nngp)

  use_par <- parallel && n_chains > 1L
  if (use_par) {
    n_workers <- max(1L, min(n_chains, parallel::detectCores(logical = FALSE) - 1L,
                             na.rm = TRUE))
    use_par <- n_workers > 1L
  }
  if (use_par) {
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(nimble))
    res <- parallel::clusterMap(
      cl, worker,
      inits_list = lapply(inits, list), seeds = as.list(seeds),
      MoreArgs = args, SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
    res <- unlist(res, recursive = FALSE)
  } else {
    # quote = TRUE stops do.call() from evaluating the model code
    res <- do.call(worker, c(list(inits_list = inits, seeds = seeds), args),
                   quote = TRUE)
  }

  coda::mcmc.list(lapply(res, function(m) {
    coda::mcmc(m, start = n_burnin + 1, thin = n_thin)
  }))
}

# Posterior processing and diagnostics ---------------------------------------

bm_diagnostics <- function(samples) {
  mat <- as.matrix(samples)
  pars <- colnames(mat)
  rhat <- rep(NA_real_, length(pars))
  if (length(samples) > 1L) {
    rhat <- tryCatch(
      coda::gelman.diag(samples, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
      error = function(e) {
        vapply(seq_along(pars), function(j) {
          tryCatch(coda::gelman.diag(samples[, j], autoburnin = FALSE)$psrf[1, 1],
                   error = function(e) NA_real_)
        }, numeric(1))
      }
    )
  }
  ess <- suppressWarnings(as.numeric(coda::effectiveSize(samples)))
  q <- t(apply(mat, 2, stats::quantile, probs = c(0.025, 0.5, 0.975)))
  data.frame(
    parameter = pars, mean = colMeans(mat), sd = apply(mat, 2, stats::sd),
    q2.5 = q[, 1], median = q[, 2], q97.5 = q[, 3],
    rhat = unname(rhat), ess = ess, row.names = NULL
  )
}

bm_convergence_warnings <- function(diag, rhat_max = 1.05, ess_min = 400) {
  if (all(is.na(diag$rhat))) return(invisible(NULL))
  bad_r <- diag[!is.na(diag$rhat) & diag$rhat > rhat_max, ]
  bad_e <- diag[!is.na(diag$ess) & diag$ess < ess_min, ]
  if (nrow(bad_r)) {
    warning(nrow(bad_r), " parameter(s) have R-hat > ", rhat_max, " (worst: ",
            bad_r$parameter[which.max(bad_r$rhat)], " = ",
            format(max(bad_r$rhat), digits = 3),
            "). Chains may not have converged; run longer.", call. = FALSE)
  }
  if (nrow(bad_e)) {
    warning(nrow(bad_e), " parameter(s) have effective sample size < ", ess_min,
            " (lowest: ", bad_e$parameter[which.min(bad_e$ess)], " = ",
            round(min(bad_e$ess)), "). Consider more iterations or thinning less.",
            call. = FALSE)
  }
  invisible(NULL)
}
