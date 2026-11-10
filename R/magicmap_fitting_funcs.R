
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
  out <- sum(sum_exp)
  return(out)
}


mstep_optim_func <- function(params, data, priorp, postp){
  
  n <- nrow(data)
  k <- ncol(postp)
  
  stopifnot(length(params)==(n+k))
  
  beta <- params[1:k]
  xpred <- params[(k+1):(n+k)]
  
  
  llik <- postp*(component_logdens(data, beta, xpred) + rep(log(priorp), each=nrow(data)))
  return(sum(llik))
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
  x_grad <- rowSums(x_grad_k)
  
  return(c(beta_grad,x_grad))
  
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
