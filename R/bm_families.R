# Response-distribution registry. Each family describes how a response enters
# the NIMBLE likelihood, so new families (e.g. Poisson) can be added here
# without changing the model builder. Besides the likelihood, a family provides
# the expected value, variance, simulation, log density and cdf used by
# prediction and the posterior predictive checks, and the name (`par2`), prior and
# initial value of its second parameter (residual sd `sigma`, or size `phi`).

bm_families <- function() {
  list(
    truncnorm = list(
      label = "truncated normal (lower = 0)",
      responses = c("biomass", "richness"),
      valid = function(y) all(y >= 0),
      valid_msg = "non-negative values",
      # Second parameter of the family (residual sd), its prior and initial value
      par2 = "sigma", par2_prior = "halft", discrete = FALSE,
      par2_init = function(scale, y) scale * stats::runif(1, 0.5, 1.5),
      # Scale on which priors are set: the response itself
      prior_scale = function(y) y,
      likelihood = function(y, mu, sigma) {
        bquote(.(y)[i] ~ dtnorm_lower(.(mu)[i], sd = .(sigma), lower = 0))
      },
      # Moments, simulation and log density of N(mu, sigma^2) truncated to
      # [0, Inf). `mu` is a draws x observations matrix; `sigma` a vector with
      # one value per draw (recycled down the columns). All work in log space
      # so they stay finite when mu is many sigmas below zero.
      expected = function(mu, sigma) {
        z <- mu / sigma
        mu + sigma * exp(stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE))
      },
      variance = function(mu, sigma) {
        z <- mu / sigma
        lam <- exp(stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE))
        sigma^2 * pmax(1 - z * lam - lam^2, 0)
      },
      simulate = function(mu, sigma) {
        z <- mu / sigma
        u <- stats::runif(length(mu))
        w <- stats::qnorm(log(u) + stats::pnorm(z, log.p = TRUE), log.p = TRUE)
        pmax(mu - sigma * w, 0)
      },
      loglik = function(y, mu, sigma) {
        stats::dnorm(y, mu, sigma, log = TRUE) - stats::pnorm(mu / sigma, log.p = TRUE)
      },
      # cdf via the survival function in log space (accurate in both tails)
      cdf = function(y, mu, sigma) {
        lsf <- stats::pnorm(y, mu, sigma, lower.tail = FALSE, log.p = TRUE) -
          stats::pnorm(mu / sigma, log.p = TRUE)
        pmin(pmax(-expm1(lsf), 0), 1)
      }
    ),
    lognormal = list(
      label = "log-normal",
      responses = "biomass",
      valid = function(y) all(y > 0),
      valid_msg = "strictly positive values",
      par2 = "sigma", par2_prior = "halft", discrete = FALSE,
      par2_init = function(scale, y) scale * stats::runif(1, 0.5, 1.5),
      prior_scale = function(y) log(y),
      likelihood = function(y, mu, sigma) {
        bquote(.(y)[i] ~ dlnorm(.(mu)[i], sdlog = .(sigma)))
      },
      expected = function(mu, sigma) exp(mu + sigma^2 / 2),
      variance = function(mu, sigma) (exp(sigma^2) - 1) * exp(2 * mu + sigma^2),
      simulate = function(mu, sigma) {
        matrix(stats::rlnorm(length(mu), mu, sigma), nrow(mu), ncol(mu))
      },
      loglik = function(y, mu, sigma) stats::dlnorm(y, mu, sigma, log = TRUE),
      cdf = function(y, mu, sigma) stats::plnorm(y, mu, sigma)
    ),
    # Poisson with a log link: y ~ Poisson(exp(mu)). No second parameter; the
    # variance equals the mean, which real richness data rarely satisfy.
    poisson = list(
      label = "Poisson (log link)",
      responses = "richness",
      valid = function(y) all(y >= 0 & y == round(y)),
      valid_msg = "non-negative integer counts",
      par2 = NULL, par2_prior = "none", discrete = TRUE,
      par2_init = function(scale, y) NULL,
      prior_scale = function(y) log1p(y),
      likelihood = function(y, mu, par2) {
        lam <- as.name(sub("^mu_", "lambda_", as.character(mu)))
        bquote({
          .(lam)[i] <- exp(.(mu)[i])
          .(y)[i] ~ dpois(.(lam)[i])
        })
      },
      expected = function(mu, par2) exp(mu),
      variance = function(mu, par2) exp(mu),
      simulate = function(mu, par2) {
        matrix(stats::rpois(length(mu), exp(mu)), nrow(mu), ncol(mu))
      },
      loglik = function(y, mu, par2) stats::dpois(y, exp(mu), log = TRUE),
      cdf = function(y, mu, par2) stats::ppois(y, exp(mu)),
      cdf_left = function(y, mu, par2) stats::ppois(y - 1, exp(mu))
    ),
    # Negative binomial with a log link: y ~ NB(mean = exp(mu), size = phi), so
    # Var(y) = exp(mu) + exp(2 mu) / phi. `mu` is the linear predictor on the log
    # scale and the family's second parameter `phi` is the size (dispersion);
    # small phi means strong overdispersion and phi -> Inf the Poisson limit.
    negbinomial = list(
      label = "negative binomial (log link)",
      responses = "richness",
      valid = function(y) all(y >= 0 & y == round(y)),
      valid_msg = "non-negative integer counts",
      par2 = "phi", par2_prior = "gamma", discrete = TRUE,
      par2_init = function(scale, y) stats::runif(1, 2, 10),
      prior_scale = function(y) log1p(y),
      likelihood = function(y, mu, phi) {
        pr <- as.name(sub("^mu_", "prob_", as.character(mu)))
        bquote({
          .(pr)[i] <- .(phi) / (.(phi) + exp(.(mu)[i]))
          .(y)[i] ~ dnbinom(prob = .(pr)[i], size = .(phi))
        })
      },
      expected = function(mu, phi) exp(mu),
      variance = function(mu, phi) {
        m <- exp(mu)
        m + m^2 / phi
      },
      simulate = function(mu, phi) {
        matrix(stats::rnbinom(length(mu), size = phi, mu = exp(mu)), nrow(mu), ncol(mu))
      },
      loglik = function(y, mu, phi) stats::dnbinom(y, size = phi, mu = exp(mu), log = TRUE),
      cdf = function(y, mu, phi) stats::pnbinom(y, size = phi, mu = exp(mu)),
      # P(Y < y), for the mid-distribution transform of discrete data
      cdf_left = function(y, mu, phi) stats::pnbinom(y - 1, size = phi, mu = exp(mu))
    )
  )
}

# Resolve the `family` argument to a named character vector (one per response)
bm_resolve_family <- function(family, responses) {
  known <- names(bm_families())
  if (is.null(family)) family <- "truncnorm"
  if (!is.character(family)) {
    stop("`family` must be a character vector.", call. = FALSE)
  }
  if (is.null(names(family))) {
    if (length(family) == 1L) {
      family <- stats::setNames(rep(family, length(responses)), responses)
    } else if (length(family) == length(responses)) {
      names(family) <- responses
    } else {
      stop("An unnamed `family` must have length 1 or one entry per response.",
           call. = FALSE)
    }
  }
  family <- family[responses]
  if (anyNA(family)) {
    stop("`family` must be named for every response: ",
         paste(responses, collapse = ", "), ".", call. = FALSE)
  }
  bad <- setdiff(family, known)
  if (length(bad)) {
    stop("Unknown family: ", paste(bad, collapse = ", "),
         ". Available: ", paste(known, collapse = ", "), ".", call. = FALSE)
  }
  fams <- bm_families()
  for (r in responses) {
    if (!r %in% fams[[family[[r]]]]$responses) {
      ok <- names(fams)[vapply(fams, function(f) r %in% f$responses, logical(1))]
      stop("Family '", family[[r]], "' is not available for response '", r,
           "'; use one of: ", paste(ok, collapse = ", "), ".", call. = FALSE)
    }
  }
  family
}
