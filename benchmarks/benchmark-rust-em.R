#!/usr/bin/env Rscript
# Time the R and rust engines on simulated effect sizes and report how closely
# they agree. Usage: Rscript benchmarks/benchmark-rust-em.R [n_variants] [nStarts]

suppressPackageStartupMessages(library(MAGICMAP))

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 200L
n_starts <- if (length(args) >= 2) as.integer(args[2]) else 10L

simulate <- function(n, slopes = c(2, 1, 0.2), props = c(0.4, 0.4, 0.2), re = 0.1, seed = 1L) {
  set.seed(seed)
  z_min <- 2
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    m <- 200000L
    true_x <- stats::rnorm(m, 0, 4)
    comp <- sample(length(slopes), m, TRUE, props)
    ex <- stats::rnorm(m)
    ey <- re * ex + sqrt(1 - re^2) * stats::rnorm(m)
    d <- data.frame(x = true_x + ex, sx = 1, y = slopes[comp] * true_x + ey, sy = 1)
    out <- rbind(out, d[abs(d$x) >= z_min, ])
  }
  out <- out[seq_len(n), ]
  out$id <- paste0("v", seq_len(n))
  out
}

# The EM fits are timed on their own (the internal magicmap_1k, chained over
# k = 1..4 as magicmap() does): a full magicmap() call also runs the SCOUTJOY
# test, which is the same for both engines and can take longer than the fits.
fit_all_k <- function(d, engine) {
  dat <- d[, c("x", "sx", "y", "sy")]
  rownames(dat) <- d$id
  prev <- NULL
  fits <- list()
  seconds <- system.time(for (k in 1:4) {
    fits[[k]] <- suppressWarnings(MAGICMAP:::magicmap_1k(
      dat, 0.1, k, nStarts = n_starts, prev = prev, engine = engine
    ))
    prev <- if (is.data.frame(fits[[k]]$slopes)) {
      list(b = fits[[k]]$slopes$b, p = fits[[k]]$component_proportions)
    }
  })[["elapsed"]]
  list(fits = fits, seconds = seconds)
}

for (seed in 1:2) {
  d <- simulate(n, seed = seed)
  r <- fit_all_k(d, "R")
  rust <- fit_all_k(d, "rust")
  slope_diff <- max(mapply(function(a, b) if (is.data.frame(a$slopes)) max(abs(a$slopes$b - b$slopes$b)) else 0,
                           r$fits, rust$fits))
  loglik_diff <- max(mapply(function(a, b) abs(a$fit_stats$logLik - b$fit_stats$logLik), r$fits, rust$fits),
                     na.rm = TRUE)
  cat(sprintf(
    "%d variants (seed %d), k = 1..4, nStarts = %d: R %.2fs, rust %.2fs (%.0fx); largest slope difference %.1e, largest logLik difference %.1e\n",
    n, seed, n_starts,
    r$seconds, rust$seconds, r$seconds / rust$seconds, slope_diff, loglik_diff
  ))
}
