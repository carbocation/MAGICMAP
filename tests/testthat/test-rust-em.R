# The rust engine is a like-for-like port of the R fit; the R code is the reference.

make_em_inputs <- function(n = 60L, slopes = c(1.5, 0.3), props = c(0.6, 0.4),
                           re = 0.2, seed = 1L) {
  set.seed(seed)
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    m <- 50000L
    true_x <- stats::rnorm(m, 0, 4)
    comp <- sample(length(slopes), m, TRUE, props)
    ex <- stats::rnorm(m)
    ey <- re * ex + sqrt(1 - re^2) * stats::rnorm(m)
    d <- data.frame(x = true_x + ex, sx = 1, y = slopes[comp] * true_x + ey, sy = 1)
    out <- rbind(out, d[abs(d$x) >= 2, ])
  }
  out <- out[seq_len(n), ]
  flip <- sign(out$x)
  out$x <- out$x * flip
  out$y <- out$y * flip
  rownames(out) <- paste0("v", seq_len(n))
  list(data = out, re = re)
}

fit_em <- function(inputs, k, engine, ...) {
  suppressWarnings(MAGICMAP:::magicmap_1k(inputs$data, inputs$re, k, engine = engine, ...))
}

expect_same_fit <- function(rust, reference, tolerance = 1e-8) {
  expect_equal(rust$fit_stats$logLik, reference$fit_stats$logLik, tolerance = tolerance)
  expect_equal(rust$fit_stats$converged, reference$fit_stats$converged)
  expect_equal(rust$fit_stats$n_starts, reference$fit_stats$n_starts)
  expect_equal(rust$fit_stats$n_starts_best, reference$fit_stats$n_starts_best)
  expect_equal(rust$slopes$b, reference$slopes$b, tolerance = tolerance)
  expect_equal(rust$slopes$se, reference$slopes$se, tolerance = 1e-6)
  expect_equal(rust$component_proportions, reference$component_proportions, tolerance = tolerance)
  expect_equal(rust$predicted_target_betas, reference$predicted_target_betas, tolerance = tolerance)
  expect_equal(as.matrix(rust$posteriors), as.matrix(reference$posteriors), tolerance = tolerance)
}

test_that("rust engine agrees with the R reference from the clustering start", {
  inputs <- make_em_inputs(seed = 11L)
  for (k in 2:4) {
    expect_same_fit(fit_em(inputs, k, "rust"), fit_em(inputs, k, "R"))
  }
})

test_that("rust engine agrees with the R reference over several starting points", {
  for (seed in c(13L, 16L)) {
    inputs <- make_em_inputs(slopes = c(2, 1, 0.2), props = c(0.4, 0.4, 0.2), seed = seed)
    previous <- fit_em(inputs, 2, "R", nStarts = 5)
    prev <- list(b = previous$slopes$b, p = previous$component_proportions)
    expect_same_fit(fit_em(inputs, 3, "rust", nStarts = 5, prev = prev),
                    fit_em(inputs, 3, "R", nStarts = 5, prev = prev))
  }
})

test_that("rust engine results do not depend on the number of threads", {
  inputs <- make_em_inputs(seed = 14L)
  one <- fit_em(inputs, 3, "rust", nStarts = 8, threads = 1)
  four <- fit_em(inputs, 3, "rust", nStarts = 8, threads = 4)
  expect_identical(one, four)
})

test_that("rust engine results do not depend on the number of threads within a start", {
  # more variants than one block of the rust engine (8192), a single start
  inputs <- make_em_inputs(n = 20000L, seed = 17L)
  one <- fit_em(inputs, 2, "rust", threads = 1)
  four <- fit_em(inputs, 2, "rust", threads = 4)
  expect_identical(one, four)
  # blocks change only the order of addition
  expect_same_fit(one, fit_em(inputs, 2, "R"))
})

test_that("rust engine is refused for estimators it does not implement", {
  inputs <- make_em_inputs(seed = 15L)
  expect_error(fit_em(inputs, 2, "rust", estimator = "numeric"), "only available")
})
