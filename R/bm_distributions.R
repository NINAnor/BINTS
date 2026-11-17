# Custom NIMBLE components evaluated in the frame that builds the model.
#
# 1. dtnorm_lower(): normal truncated below at `lower`.
#    NIMBLE's built-in T(dnorm(), lower, Inf) evaluates the normalising constant
#    as 1 - pnorm(), which underflows to 0 (log density = +Inf, accepted by
#    Metropolis samplers) when the location is far below `lower`. This version
#    computes log P(X > lower) = pnorm(., lower.tail = FALSE, log.p = TRUE), which
#    is stable in the tails.
# 2. sampler_shift: an MCMC sampler (see below).
#
# Definitions are kept as unevaluated expressions so they can be evaluated in the
# global environment (where NIMBLE looks for user-defined distributions and
# samplers), including on parallel workers that do not have this package
# installed.

bm_dist_code <- function() {
  as.call(c(as.name("{"), as.list(bm_truncnorm_code())[-1],
            as.list(bm_shift_sampler_code())[-1],
            as.list(bm_nngp_code())[-1]))
}

bm_truncnorm_code <- function() {
  quote({
    dtnorm_lower <- nimble::nimbleFunction(
      run = function(x = double(0), mean = double(0), sd = double(0),
                     lower = double(0), log = integer(0, default = 0)) {
        returnType(double(0))
        if (x < lower) {
          if (log) return(-Inf)
          return(0)
        }
        lp <- dnorm(x, mean, sd, log = TRUE) -
          pnorm(lower, mean, sd, lower.tail = FALSE, log.p = TRUE)
        if (log) return(lp)
        return(exp(lp))
      }
    )
    rtnorm_lower <- nimble::nimbleFunction(
      run = function(n = integer(0), mean = double(0), sd = double(0),
                     lower = double(0)) {
        returnType(double(0))
        # Inverse-CDF draw of Z | Z > a, using symmetry for tail accuracy
        a <- (lower - mean) / sd
        u <- runif(1, 0, 1)
        z <- -qnorm(u * pnorm(-a, 0, 1))
        x <- mean + sd * z
        if (x < lower) x <- lower
        return(x)
      }
    )
  })
}

# Custom MCMC sampler: random-walk shift of a scalar `target` (the intercept)
# together with an opposite shift of all `shiftNodes` (a spatial field), which
# leaves the likelihood unchanged. Breaks the strong posterior correlation
# between the intercept and the mean level of the field. Adaptive scale. With
# `target2` (and optionally `lamNode`), a second node is moved by
# `lambda * delta` to keep the likelihood unchanged when the field is shared.
bm_shift_sampler_code <- function() {
  quote({
    sampler_shift <- nimble::nimbleFunction(
      contains = nimble::sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        shiftNodes <- control$shiftNodes
        scale <- control$scale
        # Optional compensating node, shifted by `lambda * delta` (or `delta`),
        # for the second response of a shared latent effect
        target2 <- control$target2
        lamNode <- control$lamNode
        hasT2 <- control$hasT2
        hasLam <- control$hasLam
        calcNodes <- model$getDependencies(unique(c(target, target2, shiftNodes)))
        nAccept <- 0
        nTried <- 0
      },
      run = function() {
        lp0 <- model$getLogProb(calcNodes)
        delta <- rnorm(1, 0, scale)
        model[[target]] <<- model[[target]] + delta
        values(model, shiftNodes) <<- values(model, shiftNodes) - delta
        if (hasT2) {
          if (hasLam) {
            model[[target2]] <<- model[[target2]] + model[[lamNode]] * delta
          } else {
            model[[target2]] <<- model[[target2]] + delta
          }
        }
        lp1 <- model$calculate(calcNodes)
        jump <- (log(runif(1, 0, 1)) < lp1 - lp0)
        if (jump) {
          nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
          nAccept <<- nAccept + 1
        } else {
          nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        nTried <<- nTried + 1
        if (nTried == 100) {
          rate <- nAccept / nTried
          scale <<- scale * exp(rate - 0.44)
          nAccept <<- 0
          nTried <<- 0
        }
      },
      methods = list(reset = function() {
        nAccept <<- 0
        nTried <<- 0
      })
    )
  })
}

# NNGP latent field as ONE vector node with a custom density and ONE sampler
# (the approach of the nimbleNNGP package). A field written as 240 scalar nodes
# with dynamically indexed neighbours gives NIMBLE a huge dependency graph and
# hundreds of samplers (about 30 minutes of set-up for 240 sites); here the
# field is a single node and a compiled sampler sweeps over the sites.
#
# Density: w[1] ~ N(0, sd^2); w[l] | w[nbr(l)] ~ N(sum_j B[l,j] w[nbr[l,j]],
# sd^2 F[l]), with B and F for the length-scale grid value `idx`. `Bflat` is the
# coefficient array B[l, j, g] stored as a matrix with row (l - 1) m + j.
bm_nngp_code <- function() {
  quote({
    dnngp <- nimble::nimbleFunction(
      run = function(x = double(1), sdev = double(0), idx = integer(0),
                     Bflat = double(2), sqrtF = double(2), nbr = double(2),
                     log = integer(0, default = 0)) {
        returnType(double(0))
        n <- length(x)
        m <- dim(nbr)[2]
        lp <- dnorm(x[1], 0, sdev, log = TRUE)
        for (l in 2:n) {
          mu <- 0
          for (j in 1:m) mu <- mu + Bflat[(l - 1) * m + j, idx] * x[as.integer(nbr[l, j])]
          lp <- lp + dnorm(x[l], mu, sdev * sqrtF[l, idx], log = TRUE)
        }
        if (log) return(lp)
        return(exp(lp))
      }
    )
    rnngp <- nimble::nimbleFunction(
      run = function(n = integer(0), sdev = double(0), idx = integer(0),
                     Bflat = double(2), sqrtF = double(2), nbr = double(2)) {
        returnType(double(1))
        nn <- dim(nbr)[1]
        m <- dim(nbr)[2]
        out <- numeric(nn)
        out[1] <- rnorm(1, 0, sdev)
        for (l in 2:nn) {
          mu <- 0
          for (j in 1:m) mu <- mu + Bflat[(l - 1) * m + j, idx] * out[as.integer(nbr[l, j])]
          out[l] <- rnorm(1, mu, sdev * sqrtF[l, idx])
        }
        return(out)
      }
    )
    # Log density (and its refresh) of the nodes that depend on one site
    siteCalc_virtual <- nimble::nimbleFunctionVirtual(
      methods = list(
        calc = function() { returnType(double(0)) },
        getLP = function() { returnType(double(0)) }
      )
    )
    siteCalc <- nimble::nimbleFunction(
      contains = siteCalc_virtual,
      setup = function(model, nodes) {},
      methods = list(
        calc = function() {
          returnType(double(0))
          return(model$calculate(nodes))
        },
        getLP = function() {
          returnType(double(0))
          return(model$getLogProb(nodes))
        }
      )
    )
    # Random-walk Metropolis update of every site in turn, from its full
    # conditional: likelihood of the observations at the site, its own NNGP
    # term, and the terms of the sites that have it as a neighbour. Step sizes
    # adapt per site towards an acceptance rate of 0.44.
    sampler_nngp <- nimble::nimbleFunction(
      contains = nimble::sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        L <- control$L
        nbr <- control$nbr
        Bflat <- control$Bflat
        sqrtF <- control$sqrtF
        child <- control$child
        nChild <- control$nChild
        sdNode <- control$sdNode
        idxNode <- control$idxNode
        m <- ncol(nbr)
        stepSd <- rep(control$scale0, L)
        nAcc <- numeric(L)
        nTry <- numeric(L)
        nSweep <- 0
        calcNodes <- model$getDependencies(target)
        siteCalcs <- nimble::nimbleFunctionList(siteCalc_virtual)
        for (l in 1:L) siteCalcs[[l]] <- siteCalc(model, control$siteNodes[[l]])
      },
      run = function() {
        gi <- as.integer(model[[idxNode]])
        sdv <- model[[sdNode]]
        wv <- values(model, target)
        for (l in 1:L) {
          lp0 <- siteCalcs[[l]]$getLP() + nnTerms(l, wv, gi, sdv)
          wOld <- wv[l]
          wv[l] <- wOld + rnorm(1, 0, stepSd[l])
          values(model, target) <<- wv
          lp1 <- siteCalcs[[l]]$calc() + nnTerms(l, wv, gi, sdv)
          if (lp1 - lp0 > log(runif(1, 0, 1))) {
            nAcc[l] <<- nAcc[l] + 1
          } else {
            wv[l] <- wOld
            values(model, target) <<- wv
            dummy <- siteCalcs[[l]]$calc()
          }
          nTry[l] <<- nTry[l] + 1
        }
        dummy2 <- model$calculate(target)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        nSweep <<- nSweep + 1
        if (nSweep %% 100 == 0) {
          for (l in 1:L) {
            rate <- nAcc[l] / nTry[l]
            stepSd[l] <<- stepSd[l] * exp(rate - 0.44)
            nAcc[l] <<- 0
            nTry[l] <<- 0
          }
        }
      },
      methods = list(
        # NNGP log density terms that involve site l: its own conditional and
        # those of the sites that list it as a neighbour
        nnTerms = function(l = integer(0), wv = double(1), gi = integer(0),
                           sdv = double(0)) {
          returnType(double(0))
          lp <- 0
          if (l == 1) {
            lp <- dnorm(wv[1], 0, sdv, log = TRUE)
          } else {
            mu0 <- 0
            for (j in 1:m) mu0 <- mu0 + Bflat[(l - 1) * m + j, gi] * wv[nbr[l, j]]
            lp <- dnorm(wv[l], mu0, sdv * sqrtF[l, gi], log = TRUE)
          }
          nc <- nChild[l]
          if (nc > 0) {
            for (k in 1:nc) {
              cc <- child[l, k]
              mu1 <- 0
              for (j in 1:m) mu1 <- mu1 + Bflat[(cc - 1) * m + j, gi] * wv[nbr[cc, j]]
              lp <- lp + dnorm(wv[cc], mu1, sdv * sqrtF[cc, gi], log = TRUE)
            }
          }
          return(lp)
        },
        reset = function() {}
      )
    )
  })
}
