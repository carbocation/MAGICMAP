# Standard errors for MAGICMAP slopes and mixing proportions
#
# Uses the observed information of the mixture (incomplete-data) log-likelihood
#   l = sum_j log sum_k pi_k * phi([x_j, y_j]; [X_j, theta_k X_j], Sigma_j)
# rather than the Hessian of the EM Q function. The Q-function Hessian treats
# component memberships as known and so understates uncertainty when
# components overlap (Louis 1982).
#
# Each SNP's latent target effect X_j only enters its own term, so the
# Hessian block for X is diagonal and X can be profiled out SNP by SNP
# (Schur complement). Cost is O(n * k^2); no numerical differentiation.
#
# Parameters reported: theta_1..theta_K, pi_1..pi_{K-1}  (pi_K = 1 - sum).
#
# Arguments:
# * data: data frame with columns x, sx, y, sy, re
# * beta, priorp, xpred: fitted slopes, mixing proportions, and latent X
# * engine: "R" or "rust"; with "rust" the sums over variants are done in
#   compiled code (same terms; the R code is the reference), using `threads`
#
# Returns list with:
# * vcov: model-based covariance (inverse profiled observed information)
# * vcov_robust: sandwich covariance (I^-1 B I^-1), using per-SNP scores
# * se_slope, se_slope_robust, se_prop, se_prop_robust: SEs; se_prop has length K

magicmap_vcov <- function(data, beta, priorp, xpred, engine="R", threads=1){

  n <- nrow(data)
  K <- length(beta)
  na <- 2*K - 1          # theta_1..K, pi_1..K-1
  npar <- na + 1         # plus this SNP's X_j

  if(engine == "rust"){
    # the same per-variant terms, summed over variants in compiled code
    # (Problem::info_sums in src/rust/src/lib.rs)
    sums <- .Call(.magicmap_info_sums,
                  as.double(data$x), as.double(data$sx), as.double(data$y), as.double(data$sy),
                  as.double(data$re[1]), as.double(beta), as.double(priorp), as.double(xpred),
                  as.integer(threads))
    if(sums$counts[1] > 0){
      warning("observed information for some predicted target betas is not positive; SEs may be unreliable")
    }else if(sums$counts[2] > 0){
      # as if(any(h_xx >= 0)) does below when h_xx has NaN
      stop("missing value where TRUE/FALSE needed")
    }
    return(magicmap_vcov_from_sums(-(sums$h_aa - sums$cross), sums$b, K))
  }

  s1 <- data$sx
  s2 <- data$sy
  c12 <- data$re*s1*s2
  D <- s1^2*s2^2 - c12^2
  X <- xpred
  r1 <- data$x - X

  # per component: log(pi_k) + log density, and derivatives of the log density
  la <- matrix(NA, n, K)
  u  <- la   # d/d theta_k
  v  <- la   # d/d X
  uu <- la   # d2/d theta_k^2
  ux <- la   # d2/d theta_k dX
  vv <- la   # d2/dX^2
  for(k in 1:K){
    b <- beta[k]
    r2 <- data$y - b*X
    q <- (s2^2*r1^2 - 2*c12*r1*r2 + s1^2*r2^2)/D
    la[,k] <- log(priorp[k]) - log(2*pi) - 0.5*log(D) - 0.5*q
    u[,k]  <- -(X/D)*(c12*r1 - s1^2*r2)
    v[,k]  <- (s2^2*r1 - c12*(r2 + b*r1) + b*s1^2*r2)/D
    uu[,k] <- -s1^2*X^2/D
    ux[,k] <- -(c12*r1 - s1^2*r2 - X*c12 + b*X*s1^2)/D
    vv[,k] <- -(s2^2 - 2*b*c12 + b^2*s1^2)/D
  }

  # posteriors (stable softmax)
  m <- .row_max(la)
  P <- exp(la - m)
  P <- P/rowSums(P)

  # gradient of a_k = log pi_k + log phi_k w.r.t. (theta, pi_1..K-1, X)
  grad_k <- function(k){
    g <- matrix(0, n, npar)
    g[,k] <- u[,k]
    if(K > 1){
      if(k < K){
        g[,K+k] <- 1/priorp[k]
      }else{
        g[,(K+1):(2*K-1)] <- -1/priorp[K]
      }
    }
    g[,npar] <- v[,k]
    g
  }
  G <- lapply(1:K, grad_k)

  # per-SNP score of l_j and per-SNP Hessian (n x npar x npar)
  S <- Reduce(`+`, lapply(1:K, function(k) P[,k]*G[[k]]))
  H <- array(0, c(n, npar, npar))
  for(k in 1:K){
    Gk <- G[[k]]
    for(i in 1:npar){
      for(l in i:npar){
        H[,i,l] <- H[,i,l] + P[,k]*Gk[,i]*Gk[,l]
      }
    }
    # second derivatives of a_k
    H[,k,k]       <- H[,k,k] + P[,k]*uu[,k]
    H[,k,npar]    <- H[,k,npar] + P[,k]*ux[,k]
    H[,npar,npar] <- H[,npar,npar] + P[,k]*vv[,k]
    if(K > 1){
      if(k < K){
        H[,K+k,K+k] <- H[,K+k,K+k] - P[,k]/priorp[k]^2
      }else{
        idx <- (K+1):(2*K-1)
        for(i in idx) for(l in idx[idx >= i]) H[,i,l] <- H[,i,l] - P[,K]/priorp[K]^2
      }
    }
  }
  for(i in 1:npar){
    for(l in i:npar){
      H[,i,l] <- H[,i,l] - S[,i]*S[,l]
      H[,l,i] <- H[,i,l]
    }
  }

  # profile out each X_j
  h_xx <- H[,npar,npar]
  if(any(h_xx >= 0)){
    warning("observed information for some predicted target betas is not positive; SEs may be unreliable")
  }
  h_ax <- matrix(H[,1:na,npar], n, na)
  H_aa <- matrix(apply(H[,1:na,1:na,drop=FALSE], c(2,3), sum), na, na)
  # info = -(H_aa - sum_j h_ax h_ax' / h_xx)
  info <- -(H_aa - crossprod(h_ax, h_ax/h_xx))

  # profiled per-SNP scores for the sandwich
  S_prof <- S[,1:na,drop=FALSE] - h_ax*(S[,npar]/h_xx)
  B <- crossprod(S_prof)

  magicmap_vcov_from_sums(info, B, K)
}

# Covariances and SEs from the profiled information matrix `info` and the
# sum of outer products of the profiled per-SNP scores `B`.
magicmap_vcov_from_sums <- function(info, B, K){

  na <- 2*K - 1

  # Components that have collapsed onto the same slope make the information
  # matrix singular. Checked on the unit-free (correlation-scaled) matrix with
  # a fixed tolerance, so that whether SEs are returned for such a fit does
  # not depend on rounding error.
  info_diag <- diag(info)
  if(any(!is.finite(info)) || any(info_diag <= 0) ||
     rcond(info/sqrt(outer(info_diag, info_diag))) < 1e-10){
    stop("information matrix is singular; components may not be distinct")
  }

  vcov <- solve(info)
  vcov_robust <- vcov %*% B %*% vcov

  prop_se <- function(V){
    if(K == 1) return(0)
    Vp <- V[(K+1):na, (K+1):na, drop=FALSE]
    c(sqrt(diag(Vp)), sqrt(sum(Vp)))   # pi_K = 1 - sum(pi_1..K-1)
  }

  list(
    vcov = vcov,
    vcov_robust = vcov_robust,
    se_slope = sqrt(diag(vcov))[1:K],
    se_slope_robust = sqrt(diag(vcov_robust))[1:K],
    se_prop = prop_se(vcov),
    se_prop_robust = prop_se(vcov_robust)
  )
}
