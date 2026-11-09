
estep_posteriors <- function(data, beta, priorp, xpred){
  llik <- matrix(NA, nrow(data), length(beta))
  for(k in 1:length(beta)){
    llik[,k] <- sapply(1:nrow(data), function(a){
      log(priorp[k]) + 
        mvtnorm::dmvnorm(c(data$x[a],data$y[a]),
                         mean=c(xpred[a], beta[k]*xpred[a]),
                         sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], 
                                        data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2),
                         log=T)})
  }
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
  
  llik <- matrix(NA, nrow(data), length(beta))
  for(k in 1:length(beta)){
    llik[,k] <- sapply(1:nrow(data), function(a){
      postp[a,k] *
        (log(priorp[k]) + 
           mvtnorm::dmvnorm(c(data$x[a],data$y[a]),
                            mean=c(xpred[a], beta[k]*xpred[a]),
                            sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], 
                                           data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2),
                            log=T))})
    
  }
  return(llik)
}


ll_comp_func <- function(data, beta, priorp, xpred){
  llik <- matrix(NA, nrow(data), length(beta))
  for(k in 1:length(beta)){
    llik[,k] <- sapply(1:nrow(data), function(a){
      log(priorp[k]) + 
        mvtnorm::dmvnorm(c(data$x[a],data$y[a]),
                         mean=c(xpred[a], beta[k]*xpred[a]),
                         sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], 
                                        data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2),
                         log=T)})
  }
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
  
  
  llik <- matrix(NA, nrow(data), length(beta))
  for(k in 1:length(beta)){
    llik[,k] <- sapply(1:nrow(data), function(a){
      postp[a,k] *
        (log(priorp[k]) + 
           mvtnorm::dmvnorm(c(data$x[a],data$y[a]),
                            mean=c(xpred[a], beta[k]*xpred[a]),
                            sigma=matrix(c(data$sx[a]^2, data$re[a]*data$sx[a]*data$sy[a], 
                                           data$re[a]*data$sx[a]*data$sy[a], data$sy[a]^2),2,2),
                            log=T))})
    
  }
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
