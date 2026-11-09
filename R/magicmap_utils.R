
.get_residcor_from_call <- function(selected_model, CovIntercept, TargetXIntercept, ComparatorYIntercept){

  ###
  # extract LDSC arguments from magicmap call if possible
  ###
  if("call" %in% names(selected_model)){
    ldsc_in_call <- NULL
    if(!is.null(CovIntercept)){
      ldsc_in_call <- c(ldsc_in_call, "CovIntercept")
    }
    if(!is.null(TargetXIntercept)){
      ldsc_in_call <- c(ldsc_in_call, "TargetXIntercept")
    }
    if(!is.null(ComparatorYIntercept)){
      ldsc_in_call <- c(ldsc_in_call, "ComparatorYIntercept")
    }
    if(!is.null(ldsc_in_call)){
      warning(paste("Arguments ",paste0(ldsc_in_call,collapse=", "), "overriden by values saved in model."))
    }
    CovIntercept <- selected_model$call$CovIntercept
    TargetXIntercept <- selected_model$call$TargetXIntercept
    ComparatorYIntercept <- selected_model$call$ComparatorYIntercept
    
  }else{
    # use defaults if unspecified
    ldsc_defaulted <- NULL
    if(is.null(CovIntercept)){
      ldsc_defaulted <- c(ldsc_defaulted, "CovIntercept")
      CovIntercept <- 0
    }
    if(is.null(TargetXIntercept)){
      ldsc_defaulted <- c(ldsc_defaulted, "TargetXIntercept")
      TargetXIntercept <- 1
    }
    if(is.null(ComparatorYIntercept)){
      ldsc_defaulted <- c(ldsc_defaulted, "ComparatorYIntercept")
      ComparatorYIntercept <- 1
    }
    if(!is.null(ldsc_defaulted)){
      warning(paste(paste0(ldsc_defaulted,collapse=", "), "not specified or saved in model. Assuming no correlation."))
    }
  }
  
  re <- CovIntercept / sqrt(TargetXIntercept*ComparatorYIntercept)  # must match magicmap()

  return(re)
}




# Squared distance from each variant's observed effects to each component's
# line, in standard error units: (y - b*x)^2 / var(y - b*x). This is the
# smallest Mahalanobis distance from the point to any position on the line,
# and is chi-square with 1 df for a variant that lies on the line.
# Returns an n x k matrix.
.magicmap_line_distance <- function(x, sx, y, sy, re, slopes){
  d2 <- sapply(slopes, function(b) (y - b*x)^2/(sy^2 + b^2*sx^2 - 2*b*re*sx*sy))
  matrix(d2, nrow=length(x), ncol=length(slopes))
}

# Per-variant table of the highest-posterior component and the fit to its line.
# posteriors: data frame with x, sx, y, sy in columns 1-4 and ProbClass columns
.magicmap_assignments <- function(posteriors, slopes, re){
  pp <- as.matrix(posteriors[, grep("^ProbClass", colnames(posteriors)), drop=FALSE])
  idx <- cbind(seq_len(nrow(pp)), apply(pp, 1, which.max))
  d2 <- .magicmap_line_distance(posteriors[,1], posteriors[,2], posteriors[,3], posteriors[,4], re, slopes)[idx]
  data.frame(component = idx[,2],
             posterior = pp[idx],
             distance2 = d2,
             distance_p = stats::pchisq(d2, df=1, lower.tail=FALSE),
             row.names = rownames(posteriors))
}


#' Extract single model from MAGICMAP results object
#' 
#' @description
#' Extracts one model from a MAGICMAP results object. By default, selects the lowest-BIC model among those that converged.
#' 
#' @param model
#' output of \code{\link{magicmap}()}
#' 
#' @param which_model
#' numeric, which model to select (when multiple models were fit). Default is the lowest-BIC model among those that converged (\code{converged == TRUE}), or the lowest-BIC model overall if none converged.
#' 
#' @param hide_warnings
#' Logical, whether to skip printing warnings (highly discouraged)
#' 
#' 
#' @return
#' A list of class 'magicmap_single' containing:
#' 
#' \describe{
#'
#'   \item{\code{scoutjoy_test}}{
#'     (optional) Global test of heterogeneity from SCOUTJOY (Elliott et al., 2024). P-value is for null hypothesis that the effect sizes have a single homogeneous relationship for all variants.
#'   }
#'   
#'   \item{\code{fit_stats}}{
#'     Table of fit statistics for the selected MAGICMAP mixture model. Includes:
#'     \itemize{
#'       \item \code{logLik}: log of the complete data (marginal) likelihood
#'       \item \code{model_df}: the number of free paramters in the model, excluding nuisance parameters
#'       \item \code{BIC}: Bayesian information criterion; lower values indicate a better combination of fit and complexity
#'       \item \code{AIC}: Akaike information criterion; lower values indicate a better combination of fit and complexity
#'       \item \code{converged}: boolean indicator of whether the model fitting converged; set to NA if fitting encountered issues other than convergence (see notes field)
#'       \item \code{notes}: text description of additional considerations from fitting
#'       \item \code{iter}: number of EM iterations required for model convergence
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
#'     (optional) Data frame with one row per variant giving the mixture component with the highest posterior probability (\code{component}), that probability (\code{posterior}), and how well the variant fits that component's line: \code{distance2} is the squared distance from the observed effects to the line in standard error units, and \code{distance_p} its p-value (chi-square, 1 df). Posterior probabilities only compare the components with each other, so a variant that is far from every line can still have a posterior near 1; a small \code{distance_p} identifies such variants. Distances use the observed effects and are not adjusted for \code{TargetXPThreshold}.
#'   }
#'
#'   \item{\code{call}}{
#'     (optional) \code{\link[base]{call}} object containing all arguments used
#'   }
#'   
#'   \item{\code{profiler_metrics}}{
#'     (optional) Diagnostic metrics for each mixture component from \code{\link{get_component_diagnostics}()}
#'   }
#'   
#'   \item{\code{profiler_logic_dictionary}}{
#'     (optional) Decision criteria used to classify component homogeneity in \code{\link{get_component_diagnostics}()}
#'   }
#'
#' }
#' 
#' Unlike a full 'magicmap' object, the \code{component_proportions}, \code{slopes}, \code{predicted_target_betas}, and \code{posteriors} fields are from a single regression model, and so are not nested lists and are not named for the number of fitted components (e.g. \code{mix2components}, \code{mix3components}, etc). Mixture components in the fitted model are ordered from highest to lowest slope.
#' 
#' 
#' @export
#' 

extract_magicmap_model <- function(model, which_model=NULL, hide_warnings=FALSE){
  
  if(class(model) == "magicmap_single"){
    return(model)
  }else if(class(model) != "magicmap"){
    stop("model must be a magicmap results object")
  }
  
  if("k" %in% colnames(model$fit_stats)){
    colnames(model$fit_stats)[colnames(model$fit_stats)=="k"] <- "k_components"
  }
  

  # no selection required
  if(nrow(model$fit_stats)==1){
    
    if(is.na(model$fit_stats$converged) | !model$fit_stats$converged){
      if(!hide_warnings){
        warning("selected model failed to converge")
      }
    }
    
    k <- model$fit_stats$k_components[1]
    
  }else{
    
    if(is.null(which_model)){
      # choose among fits that converged cleanly; a collapsed or failed fit
      # (converged NA or FALSE) can have the lowest BIC but should not be chosen
      ok <- which(model$fit_stats$converged %in% TRUE)
      if(length(ok) > 0){
        which_model <- ok[which.min(model$fit_stats$BIC[ok])]
      }else{
        which_model <- which.min(model$fit_stats$BIC)
      }
    }
    
    if(!hide_warnings){
      if(is.na(model$fit_stats$converged[which_model])){
        warning("selected model had issues with fitting; check notes in fit_stats table")
      }else if(!model$fit_stats$converged[which_model]){
        warning("selected model failed to converge")
      }else if(any(is.na(model$fit_stats$converged)) | any(!model$fit_stats$converged) ){
        warning("selected model was fit successfully, but some comparison models were not")
      }
    }
    
    k <- model$fit_stats$k_components[which_model]
    
  }
  
  if(!("scoutjoy_test" %in% names(model))){
    if(!hide_warnings){
      warning("scoutjoy_test is missing; results likely old or from an internal function?")
    }
    scoutjoy_test <- data.frame(RSSobs=NA, Pvalue=NA)
  }else{
    if(model$scoutjoy_test$Pvalue > .05 && k!=1){
      if(!hide_warnings){
        warning(paste0("selected k=",k," but SCOUTJOY did not reject k=1"))
      }
    }
    scoutjoy_test <- model$scoutjoy_test
  }
  
  
    
  if(paste0("mix",k,"components") %in% names(model$posteriors)){
    
    if("call" %in% names(model)){
      out <- list(
        scoutjoy_test = scoutjoy_test,
        fit_stats = model$fit_stats[model$fit_stats$k_components==k,],
        component_proportions = model$component_proportions[[paste0("mix",k,"components")]],
        slopes = model$slopes[[paste0("mix",k,"components")]],
        predicted_target_betas = model$predicted_target_betas[[paste0("mix",k,"components")]],
        posteriors = model$posteriors[[paste0("mix",k,"components")]],
        call = model$call
      )
      
    }else{
      if(!hide_warnings){
        warning("model is missing call object; recommend rerunning magicmap if possible to save input parameters")
      }
      out <- list(
        scoutjoy_test = scoutjoy_test,
        fit_stats = model$fit_stats[model$fit_stats$k_components==k,],
        component_proportions = model$component_proportions[[paste0("mix",k,"components")]],
        slopes = model$slopes[[paste0("mix",k,"components")]],
        predicted_target_betas = model$predicted_target_betas[[paste0("mix",k,"components")]],
        posteriors = model$posteriors[[paste0("mix",k,"components")]]
      )
    }
    
  }else{
    if(length(dim(model$posteriors))==2){
      # already a single model... possibly from old magicmap_1k output that didn't differentiate class?
    
      
    }else{
      # shouldn't happen, unclear what format could end up here
      stop("model has an unexpected format; very input or retry fitting MAGICMAP?")
    }
  }
  
  if("assignments" %in% names(model)){
    out$assignments <- model$assignments[[paste0("mix",k,"components")]]
  }
  
  class(out) <- "magicmap_single"
  return(out)
  
}
