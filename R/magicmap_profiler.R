.default_magicmap_diagnostic_logic_dict <- list(
  "Insufficient N"     = "n_confident < 4",
  "Low N"              = "n_confident < 7",
  "No Spread"          = "!is.na(x_spread) && x_spread < 0.05",
  "Leveraged"          = "!is.na(lev_entropy) && lev_entropy < 0.40",
  "Somewhat Leveraged" = "!is.na(lev_entropy) && lev_entropy > 0.40 && lev_entropy < 0.60",
  "Weak Inst."         = "!is.na(mean_f) && mean_f < 50.0",
  "Fuzzy"              = "cred_ratio < 0.40",
  "Somewhat Fuzzy"     = "cred_ratio >= 0.40 && cred_ratio <= 0.60",
  "High mean D2"       = "!is.na(m1) && m1 > 2.0",
  "Overfit"            = "!is.na(ad_p) && ad_p < 0.05 && !is.na(m1) && m1 < 1.0",
  "Loose Fit"          = "!is.na(ad_p) && ad_p < 0.05 && (is.na(m1) || m1 >= 1.0)",
  "WtD2-Reject"        = "!is.na(wt_pval) && wt_pval < 0.05"
)


#' Generate Diagnostic Profiles of fitted MAGICMAP components
#'
#' @description
#' Evaluates model convergence, calculates comprehensive diagnostic metrics for each
#' component, generates a dynamic diagnostic PNG plot, and returns the updated model.
#' 
#' @param model
#' Fitted MAGICMAP model object output by \code{\link{magicmap}()}
#' 
#' @param which_model
#' numeric, which model to select (when multiple models were fit). Default is the lowest-BIC model among those that converged (\code{converged == TRUE}), or the lowest-BIC model overall if none converged.
#' 
#' @param CovIntercept 
#' For old MAGICMAP objects without saved call information, the intercept from LDSC (Bulik-Sullivan et al. 2015) genetic correlation analysis (\code{gcov_int}) used to fit MAGICMAP. 
#' 
#' @param TargetXIntercept
#' For old MAGICMAP objects without saved call information, the intercept from LDSC heritability analysis of the target (x axis) trait used to fit MAGICMAP.
#' 
#' @param ComparatorYIntercept
#' For old MAGICMAP objects without saved call information, the intercept from LDSC heritability analysis of the comparator (y axis) trait used to fit MAGICMAP.
#' 
#' @param logic_dict
#' (optional) structured list of decision criteria and flag names for defining homogeneity of each component.
#' 
#' @param verbose
#' Logical, whether to print additional logging messages
#' 
#' @return 
#' The selected MAGICMAP model, updated with profiler_metrics and profiler_logic_dictionary. See \code{\link{extract_magicmap_model}} for full format.
#'
#' @export
#' 

get_component_diagnostics <- function(model, which_model=NULL, CovIntercept=NULL, TargetXIntercept=NULL, ComparatorYIntercept=NULL, logic_dict=NULL, verbose=FALSE){
  
  if(is.null(logic_dict)){
    logic_dict <- .default_magicmap_diagnostic_logic_dict
  }
  
  # =========================================================================
  # 1. MODEL SELECTION & VALIDATION
  # =========================================================================
  
  selected_model <- extract_magicmap_model(model, which_model=which_model, hide_warnings=FALSE)
  
  n_fitted_comp <- selected_model$fit_stats$k_components
  fit_converged <- selected_model$fit_stats$converged
  
  if(verbose){
    print(sprintf("  -> Profiler evaluating best model: k=%d (Converged: %s)\n", n_fitted_comp, fit_converged))
  }
  
  best_model_params <- selected_model$slopes
  estimated_slopes <- best_model_params[, "b"]
  
  if (length(estimated_slopes) == 0) {
    stop("No estimated slopes found. Problem with model format?")
  }
  
  if(any(selected_model$component_proportions < 1e-8)){
    stop("Model includes a component with no weight. Quitting due to likely instability/poor fit.")
  }
  
  
  # =========================================================================
  # 2. EXTRACT DATA & POSTERIORS
  # =========================================================================
  # MAGICMAP posteriors inherently store the original data in columns 1-4 (x, sx, y, sy)
  full_post_df <- selected_model$posteriors
  x_vals <- full_post_df[, 1]
  sx_vals <- full_post_df[, 2]
  y_vals <- full_post_df[, 3]
  sy_vals <- full_post_df[, 4]
  
  post_data <- as.matrix(full_post_df[, -c(1:4), drop=FALSE])
  classes <- apply(post_data, MARGIN = 1, which.max)
  
  re <- .get_residcor_from_call(selected_model, CovIntercept, TargetXIntercept, ComparatorYIntercept)
  
  # raw_d2_list <- best_model_params$raw_mahalanobis_d2

    
  # =================================================================
  # 3. Compute Mahalanobis distances (D2)
  # =================================================================
  
  mahalanobis_d2 <- as.data.frame(matrix(NA, nrow=nrow(full_post_df), ncol=n_fitted_comp))
  colnames(mahalanobis_d2) <- paste0("Mahalanobis_Component",1:n_fitted_comp)
  rownames(mahalanobis_d2) <- rownames(full_post_df)
  
  

  # Pre-calculate the inverse covariance matrix elements for ALL SNPs
  # This vectorizes the math so we don't need slow loops!
  det_Sigma <- sx_vals^2 * sy_vals^2 - (re * sx_vals * sy_vals)^2
  inv_11 <- sy_vals^2 / det_Sigma
  inv_22 <- sx_vals^2 / det_Sigma
  inv_12 <- -(re * sx_vals * sy_vals) / det_Sigma
  
  for(c in 1:n_fitted_comp) {
    # 3. Calculate D^2 for ALL variants relative to this component's specific line
    resid_x <- x_vals - selected_model$predicted_target_betas
    resid_y <- y_vals - (estimated_slopes[c] * selected_model$predicted_target_betas)
    
    # The Mahalanobis quadratic form
    d2_all <- (resid_x^2 * inv_11) + (resid_y^2 * inv_22) + (2 * resid_x * resid_y * inv_12)
    
    # Save the raw vector into the list
    mahalanobis_d2[,c] <- d2_all
  }


  # Add the raw D2 list to the model
  selected_model$mahalanobis_d2 <- mahalanobis_d2

  
  # =========================================================================
  # 4. DIAGNOSTIC LOGIC & HELPER FUNCTION
  # =========================================================================

  get_cluster_diagnostics <- function(w_i, assigned_idx, raw_d2, x_vals, sx_vals) {
    n_assigned <- length(assigned_idx)
    cred_ratio <- 0; cont_cred <- 0; x_spread <- NA; lev_entropy <- NA
    m1 <- NA; m2 <- NA; ad_p <- NA; wt_d2 <- NA; wt_pval <- NA
    df_eff <- NA; mean_f <- NA
    
    n_confident <- sum(w_i > 0.95)
    eff_n <- sum(w_i)
    
    if (eff_n > 0) {
      cred_ratio <- sum(w_i[w_i > 0.95]) / eff_n
      valid_idx <- is.finite(raw_d2) & is.finite(w_i) & is.finite(x_vals) & is.finite(sx_vals)
      if (sum(valid_idx) > 0) {
        w_valid <- w_i[valid_idx]; d2_valid <- raw_d2[valid_idx]
        x_valid <- x_vals[valid_idx]; sx_valid <- sx_vals[valid_idx]
        eff_n_valid <- sum(w_valid)
        
        wt_d2 <- sum(w_valid * d2_valid)
        df_eff <- max(1e-5, eff_n_valid - 1)
        wt_pval <- pchisq(wt_d2, df = df_eff, lower.tail = FALSE)
        mean_f <- sum(w_valid * (x_valid / sx_valid)^2) / eff_n_valid
        cont_cred <- sum(w_valid^2) / eff_n_valid
      }
    }
    
    if (n_assigned > 0) {
      x_assigned <- x_vals[assigned_idx]
      sx_assigned <- sx_vals[assigned_idx]
      if (n_assigned > 1) {
        obs_variance <- var(x_assigned)
        noise_variance <- mean(sx_assigned^2)
        x_spread <- max(0, 1 - (noise_variance / obs_variance)) 
        if (n_assigned >= 4) {
          sq_dist_x <- (x_assigned)^2 
          sum_sq_x <- sum(sq_dist_x)
          if (sum_sq_x > 0) {
            safe_shares <- pmax(sq_dist_x / sum_sq_x, 1e-15) 
            lev_entropy <- (-sum(safe_shares * log(safe_shares))) / log(n_assigned) 
          } else { lev_entropy <- 0 }
        }
      } else { x_spread <- 0 }
      
      cluster_d2 <- raw_d2[assigned_idx]
      cluster_d2 <- cluster_d2[is.finite(cluster_d2)]
      if (length(cluster_d2) > 2) {
        m1 <- mean(cluster_d2); m2 <- var(cluster_d2)
        ad_p <- suppressWarnings(goftest::ad.test(cluster_d2, "pchisq", df = 1))$p.value
      }
    }
    
    if (eval(parse(text = logic_dict[["Insufficient N"]]))) {
      decision <- "Insufficient N"
    } else {
      failures <- c()
      for (rule_name in names(logic_dict)) {
        if (rule_name == "Insufficient N") next
        if (eval(parse(text = logic_dict[[rule_name]]))) failures <- c(failures, rule_name)
      }
      decision <- if (length(failures) == 0) "Currently Homogeneous" else paste(failures, collapse = " + ")
    }
    
    return(list(n_confident=n_confident, cred_ratio=cred_ratio, cont_cred=cont_cred, 
                x_spread=x_spread, lev_entropy=lev_entropy, ad_pval=ad_p, m1=m1, m2=m2, 
                wt_d2=wt_d2, wt_pval=wt_pval, df_eff=df_eff, mean_f=mean_f, diagnostic_decision=decision))
  }
  
  # =========================================================================
  # 5. COMPUTE METRICS & UPDATE MODEL OBJECT
  # =========================================================================
  diag_results <- lapply(1:n_fitted_comp, function(c) {
    get_cluster_diagnostics(w_i = post_data[, c], assigned_idx = which(classes == c), 
                            raw_d2 = selected_model$mahalanobis_d2[,c], x_vals = x_vals, sx_vals = sx_vals)
  })
  
  best_model_params <- cbind(best_model_params, dplyr::bind_rows(diag_results))
  
  selected_model$profiler_metrics <- best_model_params %>%
    dplyr::mutate(component = dplyr::row_number()) %>%
    dplyr::select(component, slope = b, n_confident, cred_ratio, cont_cred, x_spread, 
                  lev_entropy, ad_pval, m1, m2, wt_d2, wt_pval, df_eff, mean_f, diagnostic_decision)
  
  selected_model$profiler_logic_dictionary <- logic_dict
  
  return(selected_model)
  
}

