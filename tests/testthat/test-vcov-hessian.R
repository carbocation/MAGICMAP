# magicmap_vcov() derives the observed information of the mixture
# log-likelihood by hand. These tests differentiate that log-likelihood
# numerically and compare.
#
# The log-likelihood has the slopes, the mixing proportions and one latent
# target effect X_j per variant as parameters. magicmap_vcov() reports the
# covariance of (slopes, proportions) after profiling out every X_j, which is
# the inverse of
#     -(H_aa - H_ax H_xx^-1 H_xa),
# H being the Hessian over all parameters. That identity holds at any
# parameter values, so it is checked at fitted values and a little away from
# them.

# Mixture log-likelihood of all parameters, written independently of
# magicmap_vcov(): par = (slopes, pi_1..pi_{K-1}, X_1..X_n).
joint_loglik <- function(par, d, re, K) {
  slopes <- par[seq_len(K)]
  props <- if (K > 1) par[K + seq_len(K - 1)] else numeric(0)
  props <- c(props, 1 - sum(props))
  X <- par[(2 * K - 1) + seq_len(nrow(d))]
  zx <- (d$x - X) / d$sx
  dens <- 0
  for (k in seq_len(K)) {
    zy <- (d$y - slopes[k] * X) / d$sy
    quad <- (zx^2 - 2 * re * zx * zy + zy^2) / (1 - re^2)
    dens <- dens + props[k] * exp(-0.5 * quad) / (2 * pi * d$sx * d$sy * sqrt(1 - re^2))
  }
  sum(log(dens))
}

# Central-difference Hessian with one step of Richardson extrapolation.
numeric_hessian <- function(f, par, h = 1e-2, ...) {
  central <- function(h) {
    p <- length(par)
    H <- matrix(0, p, p)
    f0 <- f(par, ...)
    for (i in seq_len(p)) {
      ei <- replace(numeric(p), i, h)
      H[i, i] <- (f(par + ei, ...) - 2 * f0 + f(par - ei, ...)) / h^2
      for (j in seq_len(i - 1)) {
        ej <- replace(numeric(p), j, h)
        H[i, j] <- H[j, i] <- (f(par + ei + ej, ...) - f(par + ei - ej, ...) -
                                 f(par - ei + ej, ...) + f(par - ei - ej, ...)) / (4 * h^2)
      }
    }
    H
  }
  (4 * central(h / 2) - central(h)) / 3
}

# Profiled information for (slopes, proportions) from the numerical Hessian.
numeric_profiled_information <- function(d, re, slopes, props, X) {
  K <- length(slopes)
  na <- 2 * K - 1
  par <- c(slopes, props[seq_len(K - 1)], X)
  H <- numeric_hessian(joint_loglik, par, d = d, re = re, K = K)
  a <- seq_len(na)
  -(H[a, a, drop = FALSE] - H[a, -a, drop = FALSE] %*% solve(H[-a, -a], H[-a, a, drop = FALSE]))
}

make_vcov_inputs <- function(n, slopes, props, re, seed) {
  set.seed(seed)
  true_x <- stats::runif(n, 2, 8)
  comp <- sample(length(slopes), n, TRUE, props)
  ex <- stats::rnorm(n)
  ey <- re * ex + sqrt(1 - re^2) * stats::rnorm(n)
  sx <- stats::runif(n, 0.6, 1.4)
  sy <- stats::runif(n, 0.6, 1.4)
  data.frame(x = true_x + sx * ex, sx = sx, y = slopes[comp] * true_x + sy * ey, sy = sy,
             row.names = paste0("v", seq_len(n)))
}

# Tolerances are set by the accuracy of the numerical Hessian, whose largest
# relative error over these cases is about 1e-4 (information) and 3e-7
# (slope standard errors).
expect_information_matches <- function(d, re, slopes, props, X, standard_errors = TRUE) {
  dd <- cbind(d, re = re)
  # away from a maximum the information need not be positive definite, and
  # magicmap_vcov() then warns about the square roots it cannot take
  analytic <- suppressWarnings(MAGICMAP:::magicmap_vcov(dd, slopes, props, X))
  numeric <- numeric_profiled_information(d, re, slopes, props, X)
  expect_equal(solve(analytic$vcov), numeric, tolerance = 1e-4, ignore_attr = TRUE)
  if (standard_errors) {
    K <- length(slopes)
    expect_equal(analytic$se_slope, sqrt(diag(solve(numeric)))[seq_len(K)],
                 tolerance = 1e-5, ignore_attr = TRUE)
  }
}

cases <- list(
  list(slopes = 0.8, props = 1, re = 0),
  list(slopes = 0.8, props = 1, re = 0.4),
  list(slopes = c(1.5, 0.3), props = c(0.6, 0.4), re = 0),
  list(slopes = c(1.5, 0.3), props = c(0.6, 0.4), re = 0.4),
  list(slopes = c(1.5, 0.3), props = c(0.6, 0.4), re = -0.5),
  list(slopes = c(1.0, 0.7), props = c(0.5, 0.5), re = 0.2),      # overlapping lines
  list(slopes = c(1.5, 0.3), props = c(0.9, 0.1), re = 0.2),      # one small line
  list(slopes = c(2.0, 1.0, 0.2), props = c(0.4, 0.35, 0.25), re = 0.3)
)

test_that("analytic observed information matches the numerical Hessian at fitted values", {
  for (i in seq_along(cases)) {
    case <- cases[[i]]
    K <- length(case$slopes)
    d <- make_vcov_inputs(40L, case$slopes, case$props, case$re, seed = 100L + i)
    fit <- suppressWarnings(MAGICMAP:::magicmap_1k(d, case$re, K))
    expect_information_matches(d, case$re, fit$slopes$b, fit$component_proportions,
                               fit$predicted_target_betas)
  }
})

test_that("analytic observed information matches the numerical Hessian away from the fit", {
  for (i in seq_along(cases)) {
    case <- cases[[i]]
    K <- length(case$slopes)
    d <- make_vcov_inputs(40L, case$slopes, case$props, case$re, seed = 200L + i)
    fit <- suppressWarnings(MAGICMAP:::magicmap_1k(d, case$re, K))
    # fitted values moved off the maximum: slopes 3% steeper, the first
    # proportion 0.03 lower, every X 0.1 standard errors higher
    props <- fit$component_proportions
    if (K > 1) props <- props + c(-0.03, 0.03, numeric(K - 2))
    expect_information_matches(d, case$re, fit$slopes$b * 1.03, props,
                               fit$predicted_target_betas + 0.1 * d$sx,
                               standard_errors = FALSE)
  }
})
