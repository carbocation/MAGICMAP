
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
  out <- sum(sum_exp)
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
