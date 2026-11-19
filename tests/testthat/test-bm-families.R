test_that("truncated-normal moments match numerical integration", {
  f <- bm_families()$truncnorm
  for (m in c(-3, -0.5, 0, 1, 4)) {
    for (s in c(0.7, 2)) {
      dens <- function(y) exp(f$loglik(y, m, s))
      expect_equal(integrate(dens, 0, Inf)$value, 1, tolerance = 1e-6)
      ev <- integrate(function(y) y * dens(y), 0, Inf)$value
      va <- integrate(function(y) (y - ev)^2 * dens(y), 0, Inf)$value
      expect_equal(f$expected(matrix(m), s)[1], ev, tolerance = 1e-5)
      expect_equal(f$variance(matrix(m), s)[1], va, tolerance = 1e-4)
    }
  }
})

test_that("expected value of a truncated normal exceeds the location, and is positive", {
  f <- bm_families()$truncnorm
  mu <- matrix(c(-100, -10, -2, 0, 2, 50), 6, 1)
  e <- f$expected(mu, 3)
  expect_true(all(is.finite(e)) && all(e > 0))
  expect_true(all(e >= mu))   # equal to mu where truncation is negligible
  expect_equal(e[6], 50, tolerance = 1e-8)        # far above zero: untruncated mean
  expect_equal(e[4], 3 * sqrt(2 / pi), tolerance = 1e-8)  # mu = 0: half-normal mean
})

test_that("simulation follows the truncated distribution and stays finite in the tail", {
  f <- bm_families()$truncnorm
  set.seed(5)
  for (m in c(-2, 0.5, 3)) {
    y <- as.vector(f$simulate(matrix(m, 2e5, 1), 2))
    expect_true(all(y >= 0))
    expect_equal(mean(y), f$expected(matrix(m), 2)[1], tolerance = 0.01)
    expect_equal(stats::var(y), f$variance(matrix(m), 2)[1], tolerance = 0.03)
  }
  far <- f$simulate(matrix(-80, 1000, 1), 2)
  expect_true(all(is.finite(far)) && all(far >= 0))
  # Kolmogorov-Smirnov against the exact CDF
  m <- 1; s <- 2
  cdf <- function(y) (stats::pnorm(y, m, s) - stats::pnorm(0, m, s)) / stats::pnorm(m / s)
  y <- as.vector(f$simulate(matrix(m, 5000, 1), s))
  expect_gt(suppressWarnings(stats::ks.test(y, cdf))$p.value, 0.001)
})

test_that("log-normal family moments and simulation", {
  f <- bm_families()$lognormal
  mu <- 1; s <- 0.6
  expect_equal(f$expected(matrix(mu), s)[1],
               integrate(function(y) y * exp(f$loglik(y, mu, s)), 0, Inf)$value,
               tolerance = 1e-6)
  set.seed(2)
  y <- f$simulate(matrix(mu, 2e5, 1), s)
  expect_true(is.matrix(y) && all(y > 0))
  expect_equal(mean(y), f$expected(matrix(mu), s)[1], tolerance = 0.01)
})
