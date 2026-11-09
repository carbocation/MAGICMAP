
# Log density of each variant's observed effects under each component:
# bivariate normal for (x, y) with mean (xpred, beta_k*xpred), SDs sx and sy,
# and error correlation re. Written out in closed form so that it is evaluated
# for all variants at once. Returns an n x k matrix.
component_logdens <- function(data, beta, xpred){
  c12 <- data$re*data$sx*data$sy
  D <- data$sx^2*data$sy^2 - c12^2
  r1 <- data$x - xpred
  out <- matrix(NA, nrow(data), length(beta))
  for(k in 1:length(beta)){
    r2 <- data$y - beta[k]*xpred
    out[,k] <- -log(2*pi) - 0.5*log(D) - 0.5*(data$sy^2*r1^2 - 2*c12*r1*r2 + data$sx^2*r2^2)/D
  }
  return(out)
}

estep_posteriors <- function(data, beta, priorp, xpred){
  llik <- component_logdens(data, beta, xpred) + rep(log(priorp), each=nrow(data))
  center_ll <- llik-apply(llik,1,max)
  out <- exp(center_ll-log(rowSums(exp(center_ll))))
  # out <- signif(out, 12)
  out <- out/rowSums(out)
  return(out)
}





xpred_est_func <- function(data, beta, postp){
  numer <- matrix(NA, nrow(data), length(beta))
  denom <- matrix(NA, nrow(data), length(beta))
  covs <- data$re*data$sx*data$sy
  for(k in 1:length(beta)){
    numer[,k] <- 
      postp[,k]*(
        data$x*(data$sy^2-beta[k]*covs) + 
          data$y*(beta[k]*data$sx^2-covs)
      )
    denom[,k] <- 
      postp[,k]*(data$sy^2+beta[k]^2*data$sx^2-2*beta[k]*covs)
  }
  return(rowSums(numer)/rowSums(denom))
}

slope_est_func <- function(data, beta, xpred, postp){
  est <- rep(NA, length(beta))
  for(k in 1:length(beta)){
    w <- data$sx^2*data$sy^2*(1-data$re^2)
    numer <- (postp[,k]*xpred*(data$sx^2*data$y - data$re*data$sx*data$sy*(data$x-xpred)))/w
    denom <- (postp[,k]*data$sx^2*xpred^2)/w
    est[k] <- sum(numer)/sum(denom)
  }
  return(est)
}

slope_est_func2 <- function(data, beta, postp){
  est <- rep(NA, length(beta))
  w_ik <- matrix(NA, nrow(data), length(beta))
  ressqk <- matrix(NA, nrow(data), length(beta))
  for(k in 1:length(beta)){
    w_ik[,k] <- postp[,k]*(data$sy^2+beta[k]^2*data$sx^2-2*beta[k]*data$re*data$sx*data$sy)
    ressqk[,k] <- postp[,k]*(data$y-beta[k]*data$x)^2
  }
  wi <- rowSums(w_ik)
  ressq <- rowSums(ressqk)
  
  for(k in 1:length(beta)){
    numer <- (postp[,k]*data$x*data$y*wi-postp[,k]*data$re*data$sx*data$sy*ressq)/wi^2
    denom <- (postp[,k]*data$x^2*wi-postp[,k]*data$sx^2*ressq)/wi^2
    est[k] <- sum(numer)/sum(denom)
  }
  return(est)
} 


Q_loglik_matrix <- function(data, beta, priorp, xpred, postp){
  
  llik <- postp*(component_logdens(data, beta, xpred) + rep(log(priorp), each=nrow(data)))
  return(llik)
}


ll_comp_func <- function(data, beta, priorp, xpred){
  llik <- component_logdens(data, beta, xpred) + rep(log(priorp), each=nrow(data))
  # log-sum-exp across components, so SNPs far from every line don't underflow to -Inf
  llmax <- apply(llik,1,max)
  sum_exp <- llmax
  ok <- is.finite(llmax)
  sum_exp[ok] <- llmax[ok] + log(rowSums(exp(llik[ok,,drop=FALSE]-llmax[ok])))
  # condition on selection for the target trait (zero when there was none)
  out <- sum(sum_exp) - sum(selection_logprob(data, xpred)$logA)
  return(out)
}


mstep_optim_func <- function(params, data, priorp, postp){
  
  n <- nrow(data)
  k <- ncol(postp)
  
  stopifnot(length(params)==(n+k))
  
  beta <- params[1:k]
  xpred <- params[(k+1):(n+k)]
  
  
  llik <- postp*(component_logdens(data, beta, xpred) + rep(log(priorp), each=nrow(data)))
  return(sum(llik) - sum(selection_logprob(data, xpred)$logA))
}



mstep_gradient_func <- function(params, data, priorp, postp){
  
  n <- nrow(data)
  K <- ncol(postp)
  
  stopifnot(length(params)==(n+K))
  
  beta <- params[1:K]
  xpred <- params[(K+1):(n+K)]
  
  x_grad_k <- matrix(NA,n,K)
  beta_grad <- rep(NA,K)
  for(k in 1:K){
    beta_grad[k] <- 
      sum(-(postp[,k]/(data$sx^2*data$sy^2*(1-data$re^2)))*(
        data$re*data$sx*data$sy*(data$x-xpred)*xpred +
          data$sx^2*(beta[k]*xpred^2-data$y*xpred)
      ))
    x_grad_k[,k] <- 
      (postp[,k]/(data$sx^2*data$sy^2*(1-data$re^2)))*(
        data$sy^2*(data$x-xpred) + 
          data$re*data$sx*data$sy*(2*beta[k]*xpred-data$y-beta[k]*data$x) +
          data$sx^2*(beta[k]*data$y-beta[k]^2*xpred)
      )
  }
  x_grad <- rowSums(x_grad_k) - selection_logprob(data, xpred, derivs=TRUE)$d1
  
  return(c(beta_grad,x_grad))
  
}


# ---------------------------------------------------------------------------
# Selection on the target (x) trait
#
# When variants were included because their target-trait association passed a
# significance threshold, |x_j| >= tx_j (tx_j = z threshold * the SE used for
# selection), each variant's likelihood is conditioned on that event:
#   l_j = log sum_k pi_k * phi([x_j, y_j]; [X_j, theta_k X_j], Sigma_j) - log A_j(X_j)
#   A_j(X) = P(|x_j| >= tx_j | X) = pnorm((X - tx_j)/sx_j) + pnorm((-X - tx_j)/sx_j)
# A_j depends only on the latent target effect X_j, not on the slopes or
# mixing proportions, so the E step and the slope update given X are unchanged;
# only the log-likelihood, the update of X, and the information matrix change.
#
# data$tx holds the truncation points; when it is absent there is no selection
# and every function below returns the uncorrected result.
# ---------------------------------------------------------------------------

# log A_j(X_j) and its first two derivatives with respect to X_j
selection_logprob <- function(data, xpred, derivs=FALSE){
  n <- nrow(data)
  if(is.null(data$tx)){
    out <- list(logA=rep(0, n))
    if(derivs){
      out$d1 <- rep(0, n)
      out$d2 <- rep(0, n)
    }
    return(out)
  }
  s <- data$sx
  a <- (xpred - data$tx)/s
  b <- (-xpred - data$tx)/s
  la <- stats::pnorm(a, log.p=TRUE)
  lb <- stats::pnorm(b, log.p=TRUE)
  m <- pmax(la, lb)
  logA <- m + log(exp(la - m) + exp(lb - m))
  out <- list(logA=logA)
  if(derivs){
    # ratios phi(.)/A, computed on the log scale
    ra <- exp(stats::dnorm(a, log=TRUE) - logA)
    rb <- exp(stats::dnorm(b, log=TRUE) - logA)
    out$d1 <- (ra - rb)/s
    out$d2 <- -(a*ra + b*rb)/s^2 - out$d1^2
  }
  return(out)
}

# M-step update of the latent target effects under selection.
# For each variant the Q function in X is
#   -0.5*h*(X - X0)^2 - log A(X)
# where X0 is the update without selection (xpred_est_func) and h its
# precision. h >= 1/sx^2 and d2/dX2 log A > -1/sx^2, so this is strictly
# concave with its maximum between 0 and X0; solved by Newton steps kept
# inside that bracket.
xpred_est_func_sel <- function(data, beta, postp, tol=1e-12, maxit=100){
  x0 <- xpred_est_func(data, beta, postp)
  if(is.null(data$tx)){
    return(x0)
  }
  covs <- data$re*data$sx*data$sy
  D <- data$sx^2*data$sy^2 - covs^2
  h <- rep(0, nrow(data))
  for(k in 1:length(beta)){
    h <- h + postp[,k]*(data$sy^2 + beta[k]^2*data$sx^2 - 2*beta[k]*covs)/D
  }
  lo <- pmin(0, x0)
  hi <- pmax(0, x0)
  x <- x0
  for(i in 1:maxit){
    sel <- selection_logprob(data, x, derivs=TRUE)
    grad <- -h*(x - x0) - sel$d1
    hess <- -h - sel$d2
    # the gradient is decreasing in X, so its sign says which side the maximum is on
    lo <- ifelse(grad > 0, x, lo)
    hi <- ifelse(grad < 0, x, hi)
    x_new <- x - grad/hess
    outside <- !is.finite(x_new) | x_new < lo | x_new > hi
    x_new[outside] <- ((lo + hi)/2)[outside]
    done <- all(abs(x_new - x) <= tol*(abs(x) + data$sx))
    x <- x_new
    if(done) break
  }
  return(x)
}


# ---------------------------------------------------------------------------
# Starting points for fitting a k-component model, in addition to the
# clustering start. Fully deterministic: no random numbers are drawn.
#
# 1. From the fitted model with k-1 components (prev = list(b, p)), if given:
#    - each of its lines in turn split into two lines either side of it
#    - all of its lines plus a new line through the variant that is furthest
#      from every one of them
# 2. nStarts sets of k slopes spread over the central 90% of the per-variant
#    ratios y/x by a Halton sequence (a fixed sequence that fills a range
#    evenly; start i always has the same slopes, so a larger nStarts only
#    adds starts), with equal proportions.
#
# Returns a list of list(b, p).
# ---------------------------------------------------------------------------
magicmap_extra_starts <- function(data, k, prev=NULL, nStarts=0){
  ratios <- data$y/data$x
  ok <- is.finite(ratios)
  rng <- stats::quantile(ratios[ok], c(0.05, 0.95), names=FALSE)
  out <- list()
  
  if(!is.null(prev) && length(prev$b) == k-1 && all(is.finite(prev$b)) && all(is.finite(prev$p))){
    delta <- 0.25*diff(stats::quantile(ratios[ok], c(0.25, 0.75), names=FALSE))
    for(j in seq_along(prev$b)){
      out[[length(out)+1]] <- list(b = c(prev$b[-j], prev$b[j] - delta, prev$b[j] + delta),
                                   p = c(prev$p[-j], prev$p[j]/2, prev$p[j]/2))
    }
    d2 <- .magicmap_line_distance(data$x, data$sx, data$y, data$sy, data$re, prev$b)
    d2_nearest <- apply(d2, 1, min)
    d2_nearest[!ok] <- -Inf
    worst <- which.max(d2_nearest)
    n <- nrow(data)
    out[[length(out)+1]] <- list(b = c(prev$b, ratios[worst]),
                                 p = c(prev$p*(1 - 1/n), 1/n))
  }
  
  if(nStarts > 0){
    # first k primes as the Halton bases, one per slope
    primes <- c()
    cand <- 2
    while(length(primes) < k){
      if(all(cand %% primes != 0)) primes <- c(primes, cand)
      cand <- cand + 1
    }
    halton <- function(i, base){
      f <- 1
      r <- 0
      while(i > 0){
        f <- f/base
        r <- r + f*(i %% base)
        i <- i %/% base
      }
      r
    }
    for(i in seq_len(nStarts)){
      u <- sapply(primes, function(base) halton(i, base))
      out[[length(out)+1]] <- list(b = rng[1] + u*(rng[2] - rng[1]), p = rep(1/k, k))
    }
  }
  
  return(out)
}
