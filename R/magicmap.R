# Regression mixture model with fixed-intercept york regression components
# for a single number of components k
# 
# Arguments:
# * data: data frame with 4 columns for x, se(x), y, and se(y), in that order
# * re: scalar, correlation of sampling error between x and y, assumed to be same for all observation pairs
# * k: scalar, number of  mixture components 
# * estimator: which slope estimator to use, one of MLE, alt_MLE, or numeric
# * nStarts: number of evenly spread starting points to fit from, in addition to the clustering start and the starts built from prev; 0 fits from the clustering start only
# * prev: NULL, or list(b, p) with the slopes and proportions of the fitted model with k-1 components, used to build starting points
# * tx: NULL, or vector of truncation points for selection on the target trait (variant j was included because |x_j| >= tx_j)
# * maxIterations: maximum number of iterations for fitting; applies to both EM and M-step regression iterations
# * EMTolerance: convergence criteria for EM algorithm, in terms of relative change in loglik
# * InnerTolerance: convergence criteria for fitting slope at M step, in terms of relative change in loglik
# * verbose: whether to print loglik at each iteration; setting to >1 additional prints ll for inner slope iterations
#
# Returns:
# structured list of class 'magicmap_single', without scoutjoy_test


magicmap_1k <- function(data, re, k, estimator="MLE", maxIterations=1000, EMTolerance=1e-12, InnerTolerance=1e-12, verbose=FALSE, tx=NULL, nStarts=0, prev=NULL){
  
  # sanity checks
  if(class(data)[1] != "data.frame"){
    stop("data must be an object of class data.frame, try converting data to a data.frame \'data = as.data.frame(data)\'")
  }
  
  if(!(estimator %in% c("MLE","alt_MLE","numeric"))){
    stop("estimator must be one of \'MLE\', \'alt_MLE\', or \'numeric\'.")
  }
  
  n_starts <- 1
  n_starts_best <- 1
  
  dat <- data
  colnames(dat) <- c("x","sx","y","sy")
  dat$re <- rep(re,nrow(dat))
  
  # selection on the target trait (see selection_logprob)
  selected <- !is.null(tx)
  if(selected){
    if(estimator == "alt_MLE"){
      stop("estimator \'alt_MLE\' is not available with TargetXPThreshold; use \'MLE\' or \'numeric\'.")
    }
    dat$tx <- tx
  }
  
  ##################  
  
  # single component
  # (york regression has no selection term, so with selection k=1 is fit by the EM code below)
  if(k==1 && !selected){
    york <- suppressWarnings(yorktools::york(dat,intercept=0,method = "ML",tol=InnerTolerance,maxit=maxIterations,verbose=verbose))
    b_est <- york$b[1]
    p_est <- 1
    post_est <- matrix(1, nrow(dat), 1)
    x_est <- yorktools::pred.york(dat, model=york, intercept=0)$xpred
    ll_comp <- ll_comp_func(dat, b_est, p_est, x_est)
    em_converge <- york$converge$converged
    empty_comp <- F
    iter <- NA
    em_dec <- NA
    york_dec <- NA
    
  }else{
  
    # initial groups to get inital parameters
    ratios <- dat$y/dat$x
    ratio_var <- dat$sy^2/dat$x^2 + 
                    dat$y^2*dat$sx^2/dat$x^4 -
                    2*re*dat$y*dat$sy*dat$sx/dat$x^3
    if(k==1){
      init_clust <- rep(1, nrow(dat))
    }else{
      init_clust <- Ckmeans.1d.dp::Ckmeans.1d.dp(ratios,k=k,y=1/ratio_var)$cluster
    }
    
    # initial parameters
    b_init <- rep(NA,k)
    p_init <- rep(NA,k)
    for(i in 1:k){
      b_init[i] <- stats::lm(dat$y[init_clust==i] ~ 0 + dat$x[init_clust==i], 
                             weights=1/((dat$sx[init_clust==i]^2)+(dat$sy[init_clust==i]^2))
                             )$coefficients[1]
      p_init[i] <- mean(init_clust==i)
    }
    xpred_init <- dat$x
    
    
    b_est <- b_init
    p_est <- p_init
    x_est <- xpred_init
    ll_comp <- ll_comp_func(dat, b_est, p_est, x_est)
  
    # convergence and warning flags
    em_converge <- F
    # =========================================================================
    # BEGIN SQUAREM INTEGRATION
    # =========================================================================
    
    # 1. Define the Objective Function (Calculates the Marginal Log-Likelihood)
    em_objective <- function(theta, data, k_comps, ...) {
      b_curr <- theta[1:k_comps]
      p_curr <- theta[(k_comps+1):(2*k_comps)]
      x_curr <- theta[(2*k_comps+1):length(theta)]
      return(ll_comp_func(data, b_curr, p_curr, x_curr))
    }
    
    # 2. Define the Fixed-Point Function
    em_fixed_point <- function(theta, data, k_comps, est_method, inn_tol, max_it, ...) {
      
      # Unpack parameters
      b_curr <- theta[1:k_comps]
      p_curr <- theta[(k_comps+1):(2*k_comps)]
      x_curr <- theta[(2*k_comps+1):length(theta)]
      
      # --- E STEP ---
      post_curr <- estep_posteriors(data, b_curr, p_curr, x_curr)
      
      # --- M STEP ---
      p_new <- colMeans(post_curr)
      
      # Prevent hard zeroes which crash the likelihood
      if(any(p_new == 0)) {
        p_new <- pmax(p_new, 1e-10)
        p_new <- p_new / sum(p_new)
      }
      
      b_new <- b_curr
      x_new <- x_curr
      
      # Execute your existing M-Step estimators (The York Regression)
      if(est_method == "numeric") {
        par <- optim(c(b_curr, x_curr), 
                     fn=mstep_optim_func, gr=mstep_gradient_func,
                     data=data, priorp=pmax(p_new, 1e-16), postp=post_curr, method="BFGS",
                     control=list(fnscale=-1, maxit=max_it, reltol=inn_tol))
        b_new <- par$par[1:k_comps]
        x_new <- par$par[(k_comps+1):(nrow(data)+k_comps)]
        
      } else {
        # Internal York iterations
        ll_Q_y_old <- -Inf
        york_converge_inner <- F
        y_iter <- 0
        
        while(!york_converge_inner && y_iter < max_it) {
          y_iter <- y_iter + 1
          
          if(est_method == "MLE") {
            b_new <- slope_est_func(data, b_new, x_new, post_curr)
          } else if(est_method == "alt_MLE") {
            b_new <- slope_est_func2(data, b_new, post_curr)
          }
          
          x_new <- xpred_est_func_sel(data, b_new, post_curr)
          ll_Q_y_new <- sum(Q_loglik_matrix(data, b_new, p_new, x_new, post_curr)) - 
                          sum(selection_logprob(data, x_new)$logA)
          
          # Check inner York convergence
          if(abs((ll_Q_y_new - ll_Q_y_old) / ll_Q_y_new) < inn_tol) {
            york_converge_inner <- T
          } else {
            ll_Q_y_old <- ll_Q_y_new
          }
        }
      }
      
      # Repack and return updated parameters to SQUAREM
      return(c(b_new, p_new, x_new))
    }
    
    # 3. Fit from one starting point with the SQUAREM accelerator
    run_start <- function(b0, p0){
      sq_fit <- SQUAREM::squarem(
        par = c(b0, p0, xpred_init), 
        fixptfn = em_fixed_point, 
        objfn = em_objective,
        control = list(tol = EMTolerance, maxiter = maxIterations, trace = verbose, minimize = FALSE),
        # Arguments passed down to the fixpt and obj functions:
        data = dat, k_comps = k, est_method = estimator, inn_tol = InnerTolerance, max_it = maxIterations
      )
      b <- sq_fit$par[1:k]
      p <- sq_fit$par[(k+1):(2*k)]
      x <- sq_fit$par[(2*k+1):length(sq_fit$par)]
      list(b = b, p = p, x = x,
           # recalculated at the returned parameters: squarem's value.objfn can be from
           # an earlier evaluation when the fixed point is reached in one or two steps
           ll = ll_comp_func(dat, b, p, x),
           iter = sq_fit$fpevals, # Number of times the fixed-point function (EM cycle) was evaluated
           converge = sq_fit$convergence)
    }
    
    # 4. Fit from every starting point and keep the highest log-likelihood.
    # EM only climbs to the nearest maximum of the likelihood, and there are
    # usually several. The list of starts is fixed and ordered (no random
    # numbers; see magicmap_extra_starts), and ties go to the earlier start,
    # so the result is the same on every run. The clustering start is first.
    start_fits <- list(run_start(b_init, p_init))
    if(k > 1 && nStarts > 0){
      for(st in magicmap_extra_starts(dat, k, prev, nStarts)){
        start_fits[[length(start_fits)+1]] <- tryCatch(run_start(st$b, st$p), error = function(e) NULL)
      }
    }
    n_starts <- length(start_fits)
    
    # usable fits: finite log-likelihood and no empty component
    start_ll <- sapply(start_fits, function(f){
      if(is.null(f) || !is.finite(f$ll) || any(!is.finite(f$p)) || any(f$p < 1e-8)) -Inf else f$ll
    })
    if(any(is.finite(start_ll))){
      at_best <- start_ll >= max(start_ll) - (1e-6 + 1e-8*abs(max(start_ll)))
      chosen <- which(at_best)[1]
      n_starts_best <- sum(at_best)
    }else{
      # no usable fit: report the clustering start, as when it was the only one
      chosen <- 1
      n_starts_best <- 0
    }
    
    b_est <- start_fits[[chosen]]$b
    p_est <- start_fits[[chosen]]$p
    x_est <- start_fits[[chosen]]$x
    ll_comp <- start_fits[[chosen]]$ll
    iter <- start_fits[[chosen]]$iter
    
    # Run the E-step one final time with the optimized parameters so post_est is ready for output
    post_est <- estep_posteriors(dat, b_est, p_est, x_est)
    
    # Check if SQUAREM hit the tolerance or maxed out iterations
    # Set convergence flags expected by the rest of the script
    empty_comp <- any(p_est < 1e-8)
    em_converge <- start_fits[[chosen]]$converge
    em_dec <- NA     
    york_dec <- NA
    
    # =========================================================================
    # END SQUAREM INTEGRATION
    # =========================================================================
  }
  # model df for k-1 priors and k slopes
  df <- 2*k-1
  
  # fit stats
  BIC <- df*log(nrow(dat)) - 2*ll_comp
  AIC <- 2*df - 2*ll_comp
    
  # build output
  # order components by slope
  if(k>1){
    idx <- order(b_est, decreasing=T)
    b_est <- b_est[idx]
    p_est <- p_est[idx]
    post_est <- post_est[,idx,drop=FALSE]
  }
  
  colnames(post_est) <- paste0("ProbClass",1:k)
  rownames(post_est) <- rownames(data)
  posteriors <- cbind(data,post_est)
  
  # assigned component and fit to its line for each variant
  assignments <- .magicmap_assignments(posteriors, b_est, re)
  
  # standard errors from the observed information of the mixture likelihood
  # (see magicmap_variance.R); se is the sandwich (robust) SE, se_model the
  # model-based SE
  se_slope <- rep(NA_real_, k)
  se_model <- rep(NA_real_, k)
  if(is.finite(ll_comp) && !empty_comp){
    vc <- tryCatch(magicmap_vcov(dat, b_est, p_est, x_est),
                   error = function(e){
                     warning(paste0("Could not calculate SEs for k=", k, ": ", conditionMessage(e)))
                     NULL
                   })
    if(!is.null(vc)){
      se_slope <- vc$se_slope_robust
      se_model <- vc$se_slope
    }
  }

  parameters <- data.frame(
    b = b_est,
    se = se_slope,
    se_model = se_model
  )

  notes <- ""
  
  if(any(abs(diff(b_est,lag=1)) < 0.01)){
    warning(paste("Slopes appear collapsed to fewer than ",k,"components"))
    notes <- paste0(notes, "Slopes not well differentiated; ")
    em_converge <- NA
  }
  if(!is.na(em_dec)){
    notes <- paste0(notes, paste0("EM LL decreased ",signif(em_dec,4),"; "))
  }
  if(!is.na(york_dec)){
    notes <- paste0(notes, paste0("Inner LL decreased, max step ",signif(york_dec,4),"; "))
  }
  if(ll_comp == -Inf){
    em_converge <- NA
    notes <- paste0(notes, paste0("LL not finite; "))
  }
  
  
  out <- list(
    fit_stats = data.frame(k_components=k, 
                           logLik=ll_comp, 
                           model_df=df, 
                           BIC=BIC, 
                           AIC=AIC,
                           converged=em_converge,
                           notes=notes,
                           iter=iter,
                           n_starts=n_starts,
                           n_starts_best=n_starts_best),
    component_proportions = p_est,
    slopes = parameters,
    predicted_target_betas = x_est,
    posteriors = posteriors,
    assignments = assignments
  )
  
  if(!is.na(em_converge) && !em_converge){
    warning(paste("EM algorithm for",k,"components failed to converge in",iter,"iterations"))
  }
  
  if(empty_comp){
    warning(paste("EM algorithm for",k,"components encountered empty component"))
    em_converge <- FALSE
    
    if(notes==""){
      notes <- paste("EM halted due to empty component")
    }else{
      notes <- paste0(paste("EM halted due to empty component"),"; ",notes)
    }
    
    out <- list(
      fit_stats = data.frame(k_components=k, 
                             logLik=NA, 
                             model_df=df, 
                             BIC=NA, 
                             AIC=NA,
                             converged=F,
                             notes=notes,
                             iter=iter,
                             n_starts=n_starts,
                             n_starts_best=n_starts_best),
      component_proportions = NA,
      slopes = NA,
      predicted_target_betas = NA,
      posteriors = NA,
      assignments = NA
    )
  }
  
  class(out) <- 'magicmap_single'
  
  return(out)
  
}


#' Fit MAGICMAP model
#'
#' @description
#' Fits MAGICMAP (Mixtures Aggregating Genes Into Coordinated Modules And Pathways), a regression mixture model where each component is a york regression with the intercept fixed to zero.
#'
#' @details
#' For a full description of the model, see the vignette pdf.
#' 
#' @param data
#' Data frame with the specified column names containing genetic effect sizes and standard errors
#' 
#' @param betaTargetX 
#' String, giving the column name in \code{data} for the effect sizes (betas) on the target (x-axis) trait
#' 
#' @param sdTargetX 
#' String, giving the column name in \code{data} for the standard error (se) of the effect sizes on the target (x-axis) trait
#' 
#' @param betaComparatorY
#' String, giving the column name in \code{data} for the effect sizes (betas) on the comparator (y-axis) trait

#' @param sdComparatorY
#' String, giving the column name in \code{data} for the standard error (se) of the effect sizes on the comparator (y-axis) trait
#' 
#' @param k
#' Number of components to fit. Can be a scalar or a vector.
#' 
#' @param ids 
#' (Optional) string giving the column name in \code{data} of variant IDs to use for labeling results tables
#' 
#' @param CovIntercept 
#' Intercept from LDSC (Bulik-Sullivan et al. 2015) genetic correlation analysis (\code{gcov_int}) comparing the two traits. Default = 0, which assumes the two sets of effect sizes do not have correlated error (i.e. no sample overlap).
#' 
#' @param TargetXIntercept
#' Intercept from LDSC heritability analysis of the target (x axis) trait. Default = 1, which assumes no population stratification or other confounding.
#' 
#' @param ComparatorYIntercept
#' Intercept from LDSC heritability analysis of the comparator (y axis) trait. Default = 1, which assumes no population stratification or other confounding.
#'
#' @param forceAllK
#' Logical, whether to continue estimating the model for all values in k even if components appear to collapse at a lower value. Default false.
#' 
#' @param inflateSEs
#' Logical, whether to apply LSDC h2 intercept adjustment to the input standard errors. By default, input SEs are used as given.
#' 
#' @param TargetXPThreshold
#' (Optional) P-value threshold that every variant had to pass for the target (x axis) trait in order to be included, e.g. \code{5e-8} for genome-wide significant lead variants. When given, each variant's likelihood is conditioned on having passed that threshold (a truncated bivariate normal), which corrects the overestimation of target effects that are close to the threshold ("winner's curse") and the resulting bias of slopes towards zero. Requires that variants were selected on the target trait alone, using the betas and standard errors supplied in \code{data}; it is an error if any variant does not pass the threshold. Default NULL, which assumes variants were not selected on their target trait association. Not available with \code{estimator="alt_MLE"}. The SCOUTJOY test is not adjusted for selection.
#' 
#' @param nStarts
#' (Optional) Number of additional starting points to fit each model from. The fitting algorithm (EM) only climbs to the nearest maximum of the likelihood, and mixture likelihoods usually have several, so a single start can return a poor solution. Default 0, which fits each model from a single start (a clustering of the per-variant effect ratios). When greater than 0, each model with more than one component is fit from (1) a clustering of the per-variant effect ratios, (2) when the model with one fewer component was also fit, that solution with each line split in two and with a new line through the variant furthest from all of its lines, and (3) \code{nStarts} sets of slopes spread evenly over the range of the per-variant ratios. The fit with the highest log-likelihood is kept. No random numbers are used, so results are identical on every run. How many starts are needed depends on the data: the \code{n_starts_best} column of \code{fit_stats} gives how many starts reached the best log-likelihood, and results should be checked for stability at a larger value. Because the starts in (2) depend on the fit with one fewer component, a larger \code{nStarts} usually, but not always, gives an equal or higher log-likelihood.
#' 
#' @param estimator
#' Which estimator to use for the regression slopes. One of "MLE", "alt_MLE", or "numeric".
#' 
#' @param maxIterations
#' Maximum number of iterations for fitting each mixture model. Applies separately to EM iterations and slope fitting iterations in each M step. Default = 1000.
#' 
#' @param EMTolerance
#' Tolerance to use for EM convergence, in terms of relative change in the log likelihood. Default = 1e-12.
#' 
#' @param InnerTolerance
#' Tolerance to use for convergence of slope estimator at M step, in terms of relative change in the log likelihood. Default = 1e-12.
#' 
#' @param verbose 
#' Whether to print verbose logging. 0=no logging, 1=logging of EM iterations, 2=logging of EM and M step iterations.
#' 
#'
#' @return
#' A list of class 'magicmap' containing:
#' 
#' \describe{
#'
#'   \item{\code{scoutjoy_test}}{
#'     Global test of heterogeneity from SCOUTJOY (Elliott et al., 2024). P-value is for null hypothesis that the effect sizes have a single homogeneous relationship for all variants.
#'   }
#'   
#'   \item{\code{fit_stats}}{
#'     Table of fit statistics for the mixture model with each requested value of \code{k}. Includes:
#'     \itemize{
#'       \item \code{logLik}: log of the complete data (marginal) likelihood; conditional on selection when \code{TargetXPThreshold} is given, so not comparable to a fit without it
#'       \item \code{model_df}: the number of free paramters in the model, excluding nuisance parameters
#'       \item \code{BIC}: Bayesian information criterion; lower values indicate a better combination of fit and complexity
#'       \item \code{AIC}: Akaike information criterion; lower values indicate a better combination of fit and complexity
#'       \item \code{converged}: boolean indicator of whether the model fitting converged; set to NA if fitting encountered issues other than convergence (see notes field)
#'       \item \code{notes}: text description of additional considerations
#'       \item \code{iter}: number of EM iterations required for model convergence, for the starting point that was kept
#'       \item \code{n_starts}: number of starting points the model was fit from
#'       \item \code{n_starts_best}: how many of those reached the highest log-likelihood; when \code{nStarts} is greater than 0 and this is a small share of \code{n_starts}, a better solution may exist and a larger \code{nStarts} is worth trying
#'     }
#'   }
#'   
#'   \item{\code{component_proportions}}{
#'     Estimated prior probabilities for membership in each mixture component.
#'   }
#'   
#'   \item{\code{slopes}}{
#'     Data frame with the estimated slope (\code{b}) for each mixture component and its standard error. \code{se} is a robust (sandwich) SE and \code{se_model} is the model-based SE; both use the observed information of the mixture likelihood, so they account for uncertainty in which component each variant belongs to. \code{se} is recommended because it remains valid when slopes vary somewhat among variants within a component. SEs are NA when the fit failed.
#'   }
#'   
#'   \item{\code{predicted_target_betas}}{
#'     Estimated values of the genetic effects on the target (x axis) trait.
#'   }
#'   
#'   \item{\code{posteriors}}{
#'     Data frame containing the original input data (i.e. the X effect size, X standard error, Y effect size, and Y standard error) and the estimated posterior probability of the observations being a member of each of the \code{k} mixture components.
#'   }
#'
#'   \item{\code{assignments}}{
#'     Data frame with one row per variant giving the mixture component with the highest posterior probability (\code{component}), that probability (\code{posterior}), and how well the variant fits that component's line: \code{distance2} is the squared distance from the observed effects to the line in standard error units, and \code{distance_p} its p-value (chi-square, 1 df). Posterior probabilities only compare the components with each other, so a variant that is far from every line can still have a posterior near 1; a small \code{distance_p} identifies such variants. Distances use the observed effects and are not adjusted for \code{TargetXPThreshold}.
#'   }
#'
#'   \item{\code{call}}{
#'     \code{\link[base]{call}} object containing all arguments used 
#'   }
#'   
#' }
#' 
#' The \code{component_proportions}, \code{slopes}, \code{predicted_target_betas}, \code{posteriors}, and \code{assignments} are each structured as a list with entries named for the number of fitted components (e.g. \code{mix2components}, \code{mix3components}, etc). Within each model, mixture components are ordered from highest to lowest slope.
#' 
#' 
#' @references
#' \itemize{
#'   \item Elliott, A. \emph{et al}. Distinct and shared genetic architectures of gestational diabetes mellitus and type 2 diabetes. \emph{Nat. Genet.} 56(3), 377-382 (2024). PMCID: PMC10937370.
#'   \item Bulik-Sullivan, B. \emph{et al}. An atlas of genetic correlations across human diseases and traits. \emph{Nat. Genet.} 47(11): 1236-1241 (2015). PMCID: PMC4797329.
#' }
#'
#' @export 
#'

magicmap <- function(data, betaTargetX, sdTargetX, betaComparatorY, sdComparatorY, k, ids = NULL, CovIntercept = 0, TargetXIntercept=1, ComparatorYIntercept=1, forceAllK=FALSE, inflateSEs=FALSE, TargetXPThreshold=NULL, nStarts=0, estimator="MLE", maxIterations=1000, EMTolerance=1e-12, InnerTolerance=1e-12, verbose=FALSE){

  # sanity checks
  if((length(betaTargetX)!=1) || (length(sdTargetX)!=1) || (length(betaComparatorY)!=1) || (length(sdComparatorY)!=1)){
    stop("betaTargetX, sdTargetX, betaComparatorY, and sdComparatorY must each be a single column name")
  }
  
  if(class(data)[1] != "data.frame"){
    stop("data must be an object of class data.frame, try converting data to a data.frame \'data = as.data.frame(data)\'")
  }
  
  if(!(estimator %in% c("MLE","alt_MLE","numeric"))){
    stop("estimator must be one of \'MLE\', \'alt_MLE\', or \'numeric\'.")
  }
  
  if(!is.null(TargetXPThreshold)){
    if(!is.numeric(TargetXPThreshold) || length(TargetXPThreshold)!=1 || is.na(TargetXPThreshold) || TargetXPThreshold <= 0 || TargetXPThreshold >= 1){
      stop("TargetXPThreshold must be a single p-value threshold between 0 and 1, or NULL.")
    }
    if(estimator == "alt_MLE"){
      stop("estimator \'alt_MLE\' is not available with TargetXPThreshold; use \'MLE\' or \'numeric\'.")
    }
  }
  
  if(!is.numeric(nStarts) || length(nStarts)!=1 || is.na(nStarts) || nStarts < 0 || nStarts != round(nStarts)){
    stop("nStarts must be a single whole number, 0 or larger.")
  }
  
  k <- sort(k)
  
  # format data
  if(!is.null(ids)){
    row.names(data) <- as.character(data[,ids])
  }  
  data <- data[, c(betaTargetX, sdTargetX, betaComparatorY, sdComparatorY)]
  if(any(is.na(data))){
    na_rows <- apply(data,1,function(a) any(is.na(a)))
    warning(paste("Removing",sum(na_rows),"rows with missing values"))
    data <- data[!na_rows, ]
  }
  
  data[, c(betaComparatorY, betaTargetX)] <- data[, c(betaComparatorY, betaTargetX)] * sign(data[, betaTargetX[1]])
  
  # selection on the target trait: variants were included because |beta/se| passed
  # the z threshold, using the SEs as supplied (before any inflation below)
  tx <- NULL
  if(!is.null(TargetXPThreshold)){
    z_threshold <- stats::qnorm(TargetXPThreshold/2, lower.tail=FALSE)
    z_target <- abs(data[,betaTargetX]/data[,sdTargetX])
    if(any(z_target < z_threshold)){
      stop(paste0(sum(z_target < z_threshold), " variants do not pass TargetXPThreshold=", TargetXPThreshold,
                  " for the target trait (smallest |z| = ", signif(min(z_target),4), ", threshold |z| = ", signif(z_threshold,4),
                  "). Remove them or give the threshold that was actually used to select variants."))
    }
    tx <- z_threshold*data[,sdTargetX]
  }
  
  # handle covariance
  re <- CovIntercept / sqrt(TargetXIntercept*ComparatorYIntercept)
  
  if(inflateSEs){
    data[,sdTargetX] <- data[,sdTargetX]*sqrt(TargetXIntercept)
    data[,sdComparatorY] <- data[,sdComparatorY]*sqrt(ComparatorYIntercept)
  }
  
  scout <- SCOUTJOY::scoutjoy(betaComparatorY, betaTargetX, sdComparatorY, sdTargetX, data, CovIntercept = CovIntercept, ExposureIntercept = TargetXIntercept, OutcomeIntercept=ComparatorYIntercept, printProgress = FALSE)
  

  fits <- NULL
  slopes <- list()
  priors <- list()
  xpreds <- list()
  posts <- list()
  assigns <- list()
  stop_k_too_high <- F
  
  for(i in seq_along(k)){
    if(verbose){
      print(paste0("Testing k = ",k[i],"...")) 
    }

    # starting points for k components can be built from the fit with k-1 components
    prev <- NULL
    if(i > 1 && k[i] == k[i-1] + 1 && is.data.frame(slopes[[i-1]])){
      prev <- list(b = slopes[[i-1]]$b, p = priors[[i-1]])
    }
    
    mod <- magicmap_1k(data=data, 
                       re=re, 
                       k=k[i], 
                       estimator=estimator, 
                       maxIterations=maxIterations, 
                       EMTolerance=EMTolerance, 
                       InnerTolerance=InnerTolerance, 
                       verbose=verbose,
                       tx=tx,
                       nStarts=nStarts,
                       prev=prev)
    
    fits <- rbind(fits, mod$fit_stats)
    priors[[i]] <- mod$component_proportions
    slopes[[i]] <- mod$slopes
    xpreds[[i]] <- mod$predicted_target_betas
    posts[[i]] <- mod$posteriors
    assigns[[i]] <- mod$assignments
    
    if(any(grepl("Slopes not well differentiated",mod$fit_stats$notes))){
      if(forceAllK){
        warning(paste0("Components collapsed at k=",k[i]))
      }else{
        stop_k_too_high <- T
        if(i < length(k)){
          k <- k[1:i]
          warning(paste0("Components collapsed at k=",k[i],"; skipping remaining k values. Set forceAllK=TRUE to avoid this behavior."))
        }
        break
      }
    }
    
        
  }

  names(priors) <- paste0("mix",k,"components")
  names(slopes) <- paste0("mix",k,"components")
  names(xpreds) <- paste0("mix",k,"components")
  names(posts) <- paste0("mix",k,"components")
  names(assigns) <- paste0("mix",k,"components")

  # record call
  call_obj <- rlang::call_match(defaults=TRUE)
  fn <- call_obj[[1]]
  data_name <- call_obj$data
  args <- rlang::call_args(call_obj)
  # evaluate arguments where magicmap() was called from; eval.parent inside
  # sapply() would look in the global environment instead
  caller_env <- parent.frame()
  arg_vals <- lapply(args[-1], eval, envir = caller_env)
  clean_call <- as.call(c(fn, data_name, arg_vals))

  out <- list(
    scoutjoy_test=scout$`Global Test`,
    fit_stats = fits,
    component_proportions = priors,
    slopes = slopes,
    predicted_target_betas = xpreds,
    posteriors = posts,
    assignments = assigns,
    call = clean_call
  )
  
  class(out) <- 'magicmap'
  
  if(scout$`Global Test`$Pvalue > .05){
    warning("!! Scoutjoy global test is not significant, cannot reject 1 component solution !!")
    out$fit_stats$notes <- paste("SCOUTJOY does not reject k=1;",out$fit_stats$notes)
  }
  
  return(out)
  
}
