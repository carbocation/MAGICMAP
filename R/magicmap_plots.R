.magicmap_default_colors <- function(k){
  
  if(k > 19){
    stop(paste("Default colors are only available for k <= 19. Please specify colors."))
  }
  
  if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
    colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d",
                "#e41a1c","#377eb8","#984ea3","#ffff33","#a65628",
                "#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494")
  }else{
    colors <- c(RColorBrewer::brewer.pal(7,"Dark2"), RColorBrewer::brewer.pal(7,"Set1")[c(1,2,4,6,7)], RColorBrewer::brewer.pal(7,"Set2"))
  }
  
  colors <- c(colors[1:k],"gray70")
}



#' Basic Scatter Plot of MAGICMAP results
#' 
#' @description
#' Makes a scatter plot showing the regression mixture components fitted to the two sets effect size estimates by a MAGICMAP model.
#' 
#' 
#' @param model
#' output of \code{\link{magicmap}()}
#' 
#' @param which_model
#' numeric, which model to plot (when multiple models were fit). Default is the model with lowest BIC.
#' 
#' @param class_thresh
#' threshold of posterior probability to use for assigning variants to components
#' 
#' @param diagnostics
#' Which homogeneity diagnostics to annotate. Can include "summary" (legend with categorical decisions) and/or "details" (metrics annotated on regression line).
#' 
#' @param outer_title
#' (optional) string for a title added to the outer margin. If specified, SCOUTJOY global test results are also appended to the title. 
#' 
#' @param label_comp
#' Numeric, which mixture component(s) to label points from. Unclassified points are treated as the k+1 component.
#' 
#' @param colors
#' set of k+1 colors to use for plotting. Last color is used for unclassified variants.
#' 
#' @param se_bars
#' logical, whether to include standard error bars
#' 
#' @param se_color
#' color(s) to use for standard errors, either a single color or a vector of k+1 colors
#' 
#' @param legend
#' logical, whether to include a legend
#' 
#' @param comp_names
#' list of names to use for mixture components (only affects legend)
#' 
#' @param se_length
#' length for caps on SE intervals; passed to \code{\link[graphics]{arrows}()}
#' 
#' @param conf_region
#' Logical, whether to include 95\% confidence bands for slopes based on their SEs 
#' 
#' @param assigned_cex
#' point size scaling factor for variants assigned to a component
#' 
#' @param unassigned_cex
#' point size scaling factor for variants that are not assigned to a component
#' 
#' @param label_cex
#' scaling factor for text labels (if applicable)
#' 
#' @param CovIntercept 
#' For adding profiling to old MAGICMAP objects without saved call information, the intercept from LDSC (Bulik-Sullivan et al. 2015) genetic correlation analysis (\code{gcov_int}) used to fit MAGICMAP. 
#' 
#' @param TargetXIntercept
#' For adding profiling to old MAGICMAP objects without saved call information, the intercept from LDSC heritability analysis of the target (x axis) trait used to fit MAGICMAP.
#' 
#' @param ComparatorYIntercept
#' For adding profiling to old MAGICMAP objects without saved call information, the intercept from LDSC heritability analysis of the comparator (y axis) trait used to fit MAGICMAP.
#'
#' @param hide_warnings
#' Logical, whether to skip printing warnings (highly discouraged)
#'
#' @param ...
#' additional parameters to be passed directly to \code{\link[base]{plot}()}
#'
#' @export plot.magicmap
#' @export
#'

plot.magicmap <- function(model, which_model=NULL, class_thresh=0.95, diagnostics=c("summary"), outer_title=NULL, label_comp=NULL, colors=NULL, se_bars=TRUE, se_color="gray70", legend=TRUE, comp_names=NULL, se_length=0.025, conf_region=TRUE, assigned_cex=1, unassigned_cex=1, label_cex=0.8, CovIntercept=NULL, TargetXIntercept=NULL, ComparatorYIntercept=NULL, hide_warnings=FALSE, ...){  
  args <- list(...)
  
  selected_model <- extract_magicmap_model(model, which_model=which_model, hide_warnings=hide_warnings)
  
  x <- selected_model$posteriors[,1]
  y <- selected_model$posteriors[,3]
  sx <- selected_model$posteriors[,2]
  sy <- selected_model$posteriors[,4]
  
  if(!("xlab" %in% names(args))){
    args$xlab <- colnames(selected_model$posteriors)[1]
  }
  
  if(!("ylab" %in% names(args))){
    args$ylab <- colnames(selected_model$posteriors)[3]
  }
  

  k <- selected_model$fit_stats$k_components[1]
  
  slopes <- selected_model$slopes[,"b"]
  slopese <- selected_model$slopes[,"se"]
  
  
  # note: currently unused
  has_noise <- "ProbNoise" %in% colnames(selected_model$posteriors)
  if (has_noise) {
    noise_idx <- k + 1
    unassigned_idx <- k + 2
  } else {
    unassigned_idx <- k + 1
  }
  
  if(k==1 && !has_noise){
    classes=rep(1, nrow(selected_model$posteriors))
  }else{
    classes <- apply(selected_model$posteriors[,-c(1:4)], MARGIN = 1, which.max)
    maxprob <- apply(selected_model$posteriors[,-c(1:4)], MARGIN = 1, max)
    classes[maxprob < class_thresh] <- unassigned_idx
  }
  
  unassigned_label <- sprintf("Unassigned (N=%d)", sum(classes == unassigned_idx))
  if(has_noise){
    noise_label <- sprintf("Unstructured (N=%d)", sum(classes == noise_idx))
  }else{
    noise_label <- NULL
  }
  if(is.null(comp_names)){
    comp_labels <- c(sapply(1:k, function(a) sprintf("Component %s (N=%d)", a, sum(classes == a))))
    comp_labels <- c(comp_labels, unassigned_label, noise_label)
  }else if(length(comp_names)==k){
    comp_labels <- c(sapply(1:k, function(a) sprintf("%s (N=%d)", comp_names[a], sum(classes == a))))
    comp_labels <- c(comp_labels, unassigned_label, noise_label)
  }else{
    comp_labels <- c(sapply(1:length(comp_names), function(a) sprintf("%s (N=%d)", comp_names[a], sum(classes == a))))
  }
  

  if(any(!(diagnostics %in% c("summary","details")))){
    stop("Valid options for diagnostics are 'summary' and/or 'details' or NULL.")
  }else if(!is.null(diagnostics)){
    if(!("profiler_metrics" %in% names(selected_model))){
      selected_model <- get_component_diagnostics(selected_model, CovIntercept=CovIntercept, TargetXIntercept=TargetXIntercept, ComparatorYIntercept=ComparatorYIntercept)
    }
  }

  if("details" %in% diagnostics){
    slopes_annot <- selected_model$profiler_metrics %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        annotation_font_size = 0.4,
        ad_safe = ifelse(is.na(ad_pval), "NA", sprintf("%.3f", ad_pval)),
        lev_safe = ifelse(is.na(lev_entropy), "NA", sprintf("%.2f",lev_entropy)),
        wt_p_safe = ifelse(is.na(wt_pval), "NA", sprintf("%.3f", wt_pval)),
        label_text = sprintf("N_conf: %d | IGX2: %.2f | LevEnt: %s\nCred: %.2f | wtF-stat: %.1f |\nAD p: %s | m1: %.2f\nwtD2: %.2f (df=%.1f,p=%s)",
                             n_confident, x_spread, lev_safe, cred_ratio, mean_f, ad_safe, 
                             m1, wt_d2, df_eff, wt_p_safe)
      ) %>% dplyr::ungroup()
    # selected_model$profiler_metrics_formatted <- slopes_annot
  }
  
  if("summary" %in% diagnostics){
    legend <- FALSE
    wrapped_decisions <- sapply(selected_model$profiler_metrics$diagnostic_decision, function(text) {
      paste(strwrap(text, width = 25), collapse = "\n  ")
    })
    custom_legend <- sprintf("%s:\n  %s", comp_labels[seq_len(k)], wrapped_decisions)

    if (has_noise) {
      custom_legend <- c(custom_legend, noise_label, unassigned_label)
      leg_lty <- c(rep(1, k), NA, NA)  
      leg_lwd <- c(rep(3, k), NA, NA)  
      leg_pch <- c(rep(NA, k), 15, 15) 
    } else {
      custom_legend <- c(custom_legend, unassigned_label)
      leg_lty <- c(rep(1, k), NA)  
      leg_lwd <- c(rep(3, k), NA)  
      leg_pch <- c(rep(NA, k), 15) 
    }
    
    total_legend_items <- length(custom_legend)
    leg_y_gap <- ifelse(total_legend_items > 6, 1.4, 2.5)
    leg_txt_size <- ifelse(total_legend_items > 6, 0.5, 0.7)
    leg_font <- ifelse(total_legend_items > 6, 1, 2) 
  }
  
  
  if(is.null(colors)){
    if(has_noise){
      colors <- c(.magicmap_default_colors(k), "#666666")
    }else{
      colors <- .magicmap_default_colors(k)
    }
  }
  
  if(!(length(se_color) %in% c(1, k+1))){
    stop("se_color must either a single value or k+1 values")
  }
  
  if(!("xlim" %in% names(args))){
    args$xlim <- c(min(0,min(x-sx)), 1.25*max(x+sx))
  }
  
  if(!("xaxs" %in% names(args))){
    args$xaxs <- "i"
  }
  
  if(!("ylim" %in% names(args))){
    args$ylim <- c(min(0,min(y-sy)), 1.25*max(y+sy))
  }
  
  if(!("main" %in% names(args))){
    args$main <- ""
  }
  
  if(!("pch" %in% names(args))){
    args$pch <- 20
  }
  
  if("cex" %in% names(args)){
    if(!hide_warnings){
      warning("cex is not used; see assigned_cex, unassigned_cex, and label_cex or more specific default parameters (e.g. cex.axis, cex.lab)")
    }
  }
  
  if(!("cex.axis" %in% names(args))){
    args$cex.axis <- 0.8
  }
  
  if(!("col" %in% names(args))){
    args$col <- colors[classes]
  }
  
  if(!("lwd" %in% names(args))){
    args$lwd <- 2
  }
  
  if(!("las" %in% names(args))){
    args$las <- 1
  }
  
  if(!("tcl" %in% names(args))){
    args$tcl <- -0.25
  }
  
  if(!("mar" %in% names(args))){
    mar <- c(3.5, 
             3.5, 
             ifelse(nchar(args$main)>0, 3.5, 0.25), 
             ifelse("summary" %in% diagnostics, 0.25, 0.25)
    )
  }else{
    mar <- args$mar
  }
  
  if(!("mgp" %in% names(args))){
    mgp <- c(2.5,0.5,0)
  }else{
    mgp <- args$mgp
  }
  
  mar_bak <- par()$mar
  mgp_bak <- par()$mgp
  
  if(!is.null(outer_title)){
    oma_bak=par()$oma
    par(oma = c(0, 0, 5, 0), mar=mar, mgp=mgp)
  }else{
    par(mar=mar, mgp=mgp)
  }
  
  args$x <- 0
  args$y <- 0
  plotcol <- args$col
  args$col <- rgb(0,0,0,0)
  args$bty='l'
  
  if("summary" %in% diagnostics){
    nf <- layout(matrix(c(1,2), ncol = 2), widths = c(2/3,1/3), heights = 1)
  }
  do.call(plot, args)
  args$col <- plotcol
  
  if(conf_region){
    for(i in 1:k){
      if(!is.na(slopese[i])) {
        
        if(min(x-sx)<0){
          polygon(x=c(
            0,
            min(x-sx), 
            min(x-sx),
            0,
            max(x+sx), 
            max(x+sx)
          ),
          y=c(
            0,
            min(x-sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i]), 
            min(x-sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            0,
            max(x+sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            max(x+sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i])
          ),
          border=FALSE,
          col=adjustcolor(colors[i], alpha.f=0.3))
        }else{
          polygon(x=c(
            min(x-sx), 
            min(x-sx),
            max(x+sx), 
            max(x+sx)
          ),
          y=c(
            min(x-sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i]), 
            min(x-sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            max(x+sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            max(x+sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i])
          ),
          border=FALSE,
          col=adjustcolor(colors[i], alpha.f=0.3))
        }
        
      }
    }
  }
  
  abline(h=0,col="grey20",lwd=0.5)
  abline(v=0,col="grey20",lwd=0.5)
  abline(0,1,col="grey80",lty=2)
  abline(0,-1,col="grey80",lty=2)
  
  
  if(se_bars){
    if(length(se_color)==1){
      arrows(x0=x, x1=x, y0=y-sy, y1=y+sy, length=se_length, angle=90, code=3, col=se_color)
      arrows(x0=x-sx, x1=x+sx, y0=y, y1=y, length=se_length, angle=90, code=3, col=se_color)
    }else{
      arrows(x0=x, x1=x, y0=y-sy, y1=y+sy, length=se_length, angle=90, code=3, col=se_color[classes])
      arrows(x0=x-sx, x1=x+sx, y0=y, y1=y, length=se_length, angle=90, code=3, col=se_color[classes])
    }
  }
  
  points(x[classes==(k+1)], y[classes==(k+1)],
         pch=args$pch,
         col=args$col[classes==(k+1)],
         cex=unassigned_cex)
  
  for(i in 1:k){
    points(x[classes==i], y[classes==i],
           pch=args$pch,
           col=args$col[classes==i],
           cex=assigned_cex)
  }
  
  for(i in 1:k){
    abline(0, slopes[i], col=colors[i], lwd=args$lwd)
  }
  
  
  if(legend){

    legend("topright", 
           fill=colors,
           border=colors,
           # lwd=rep(lwd,k+1), 
           cex = label_cex, 
           bg="white",
           box.col="white", 
           legend = comp_labels)
  }
  

  if(!is.null(label_comp)){
    
    df_post <- selected_model$posteriors

    
    for(c in label_comp){
      text(x=df_post[classes==c, 1], 
           y=df_post[classes==c, 3], 
           labels=rownames(df_post)[classes==c], 
           cex=label_cex, 
           adj=c(-.1,-.25))
    }
  }
  
  if("details" %in% diagnostics){
    par(xpd = NA) 
    usr <- par("usr"); xmin=usr[1]; xmax=usr[2]; ymin=usr[3]; ymax=usr[4]
    
    mx <- (xmax - xmin) * 0.02; my <- (ymax - ymin) * 0.02
    
    for (j in seq_len(nrow(slopes_annot))) {
      cand_x <- seq(max(0, xmin), xmax, length.out=1000); cand_y <- slopes_annot$slope[j] * cand_x
      valid <- which(cand_y >= ymin & cand_y <= ymax)
      if(length(valid)>0){ idx <- valid[max(1, round(length(valid)*0.85))]; ax=cand_x[idx]; ay=cand_y[idx]
      } else { ax=(xmax+max(0,xmin))/2; ay=(ymax+ymin)/2 }
      
      adj_x <- 1; adj_y <- ifelse(ay > (ymin + ymax)/2, 1, 0)
      label_cex <- slopes_annot$annotation_font_size[j]
      tw <- strwidth(slopes_annot$label_text[j], cex=label_cex, font=2)
      th <- strheight(slopes_annot$label_text[j], cex=label_cex, font=2)
      
      bx1 <- ax - tw; bx2 <- ax; by1 <- ifelse(adj_y == 1, ay - th, ay); by2 <- ifelse(adj_y == 1, ay, ay + th)
      if (bx2 > xmax - mx) ax <- ax - (bx2 - (xmax - mx))
      if (bx1 < xmin + mx) ax <- ax + ((xmin + mx) - bx1)
      if (by2 > ymax - my) ay <- ay - (by2 - (ymax - my))
      if (by1 < ymin + my) ay <- ay + ((ymin + my) - by1)
      
      text(ax, ay, labels=slopes_annot$label_text[j], adj=c(adj_x, adj_y), cex=label_cex, font=2, col=colors[j])
    }
    par(xpd = FALSE)
  }
  
  if("summary" %in% diagnostics){
     
    par(mar=c(0,0,mar[3],mar[4]))
    plot(0,0,
         col=rgb(0,0,0,0),
         bty='n',
         xaxt='n',
         yaxt='n',
         xlab="",
         ylab="",
        
         )
    par(xpd = NA)
    legend("left", 
           # inset = c(-.35, 0), 
           legend = custom_legend, 
           col = colors, 
           lty = leg_lty, 
           lwd = leg_lwd, 
           pch = leg_pch, 
           pt.cex = leg_txt_size * 2.2, 
           text.col = colors, 
           bty = "n", 
           cex = leg_txt_size, 
           text.font = leg_font, 
           y.intersp = leg_y_gap,
           bg="white",
           box.col="white"
    )
    
    par(mar=mar)
    nf <- layout(matrix(1, ncol = 1), widths = 1, heights = 1)
    par(xpd = FALSE)
  }
  
  
  if(!is.null(outer_title)){
    sj_pval <- model$scoutjoy_test$Pvalue
    sj_reject <- ifelse(sj_pval < 0.05, "TRUE", "FALSE")
    title_line1 <- sprintf("%s | SCOUTJOY Reject: %s (p=%.2g)", outer_title, sj_reject, sj_pval)
    
    base_cex <- 1.2
    max_width_allowed <- 0.7 
    current_width <- strwidth(title_line1, units="figure", cex=base_cex)
    fit_cex <- if(current_width > max_width_allowed) base_cex * (max_width_allowed / current_width) else base_cex
    fit_cex <- max(fit_cex, 0.6)
    
    mtext(title_line1, outer = TRUE, side = 3, line = 1, font = 2, cex = fit_cex, adj = 0.5) 
    par(mar=mar_bak, mgp=mgp_bak, oma=oma_bak)
  }else{
    par(mar=mar_bak, mgp=mgp_bak)
  }

}

#' @rdname plot.magicmap
#' @export
plot.magicmap_single <- plot.magicmap



#' MAGICMAP results figure showing distribution of ratio of effect sizes
#' 
#' @description
#' Alternative plot of MAGICMAP results based on the estimated distribution of ratios between the two input effect sizes. By default includes a histogram of the estimated per-component distributions, a forest plot of the estimated ratio for each variant, and a scatter plot showing the fitted regression.
#' 
#' @param model
#' output of \code{\link{magicmap}}()
#' 
#' @param target_name
#' string, name for labeling the target (x axis) trait
#' 
#' @param comparator_name
#' string, name for labeling the comparator (y axis) trait
#' 
#' @param panels
#' which panels to include in the figure. Options are: "all","hist","histforest" (histogram+forest plot),"forest", or "foresttable" (forest plot+table)
#' 
#' @param output_forest_data
#' whether to output the "raw" plotted data, the summarized values for the "forest" plot, or "none" (default) 
#' 
#' @param which_model
#' numeric, which model to plot (when multiple models were fit). Default is the model with lowest BIC.
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
#' @param class_thresh
#' threshold of posterior probability to use for assigning variants to components
#' 
#' @param forest_drop_unassigned
#' omit unassigned points from forest plot (if applicable)
#' 
#' @param conf_region
#' Logical, whether to include 95\% confidence bands for slopes based on their SEs 
#' 
#' @param colors
#' set of k+1 colors to use for plotting. Last color is used for unclassified variants.
#' 
#' @param cex_forest_labels
#' size scaling for text in forest plot/table
#' 
#' @param layout_widths
#' relative widths for the columns of the plot when plotting "all" panels or "foresttable" (forest plot+table). Default = c(3,2).
#' 
#' @param layout_heights
#' relative heights for the rows of the plot when plotting "all" panels or "histforest" (histogram+forest plot). Default = c(1,3).
#' 
#' @param n_mc
#' number of monte carlo replicates to simulate to estimate the distribution of effect size ratios
#' 
#' @param seed
#' seed for monte carlo estimation of distribution of effect size ratios
#'
#' @param ...
#' additional parameters (currently ignored)
#'
#'
#' @export
#'

plot_ratio_figure <- function(model, target_name=NULL, comparator_name=NULL, panels="all", output_forest_data="none", which_model=NULL, CovIntercept=NULL, TargetXIntercept=NULL, ComparatorYIntercept=NULL, class_thresh=0.95, forest_drop_unassigned=FALSE, conf_region=TRUE, colors=NULL, cex_forest_labels=0.5, layout_widths=NULL, layout_heights=NULL, n_mc=20000, seed=NULL, ...){  
  
  args <- list(...)
  
  if(class(model) != "magicmap" && class(model) != "magicmap_single"){
    stop("model must be a magicmap results object")
  }
  
  ##########
  # Setup
  ##########
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  num_breaks = 501
  
  panels_opts <- c("all","hist","histforest","forest","foresttable")
  if(!(panels %in% panels_opts)){
    stop(paste("panels must be one of:",paste(panels_opts,collapse = ",")))
  }
  
  if(!(output_forest_data %in% c("raw","forest","none"))){
    stop(paste("output_forest_data must be one of:",paste(c("raw","forest","none"),collapse = ",")))
  }
  
  
  # for old yorkmix output, fake output format
  if(class(model) != "magicmap" && class(model) != "magicmap_single"){
    warning("Does not appear to be magicmap output, will try to convert")
    
    model$slopes <- model$parameters
    for(i in 1:nrow(model$fit_stats)){
      model$slopes[[i]] <- cbind(model$slopes[[i]],se=NA)
    }
    model$scoutjoy_test <- data.frame(Pvalue=0,RSS=NA)
    class(model) <- 'magicmap'
    
  }
  
  ###
  # get info about model to plot
  ###
  
  selected_model <- extract_magicmap_model(model, which_model=which_model, hide_warnings=FALSE)

  if(is.null(target_name)){
    target_name <- colnames(selected_model$posteriors)[1]
  }
  if(is.null(comparator_name)){
    comparator_name <- colnames(selected_model$posteriors)[3]
  }

  k <- selected_model$fit_stats$k_components[1]
  tab <- selected_model$posteriors
  
  rr <- .get_residcor_from_call(selected_model, CovIntercept, TargetXIntercept, ComparatorYIntercept)

  
  ###
  # colors
  ###
  if(is.null(colors)){
    colors <- .magicmap_default_colors(k)
  }

  
  
  ##########
  # Estimate distribution of ratios
  ##########
  
  # simulate distribution of ratio for each variant
  mc <- data.frame(id=rep(rownames(tab), each=n_mc),
                   ratio=NA)
  
  for(i in 1:nrow(tab)){
    sim <- mvtnorm::rmvnorm(n_mc,
                            mean = as.numeric(tab[i,c(1,3)]),
                            sigma = matrix(c(tab[i,2]^2,
                                             rr*tab[i,2]*tab[i,4],
                                             rr*tab[i,2]*tab[i,4],
                                             tab[i,4]^2),2,2))
    mc$ratio[mc$id==rownames(tab)[i]] <- sim[,2]/sim[,1]
  }
  
  
  # summarize CI, group for each snp
  forest <- data.frame(
    id=rownames(tab),
    rat=tab[,3]/tab[,1],
    lb=NA,
    ub=NA,
    cat=NA,
    prob=NA
  )
  
  for(i in 1:nrow(tab)){
    cat <- which(tab[i,5:ncol(tab)] > class_thresh)
    if(length(cat)==0){
      cat <- k+1
      forest$prob[i] <- max(tab[i,5:ncol(tab)])
    }else{
      cat <- which.max(tab[i,5:ncol(tab)])
      forest$prob[i] <- tab[i,c(5:ncol(tab))[cat]]
    }
    forest$cat[i] <- cat
    
    forest$lb[i] <- quantile(mc$ratio[mc$id==rownames(tab)[i]], .025)
    forest$ub[i] <- quantile(mc$ratio[mc$id==rownames(tab)[i]], .975)
  }
  
  cols <- colors[1:(k+1)]
  forest$col <- colors[forest$cat]

  
  
  ##########
  # Begin plotting
  ##########

  ###
  # setup layout
  ###
  par_backup <- par(no.readonly = TRUE)
  
  if(panels=="all"){
    if(is.null(layout_widths)){
      widths <- c(3,2)
    }else{
      if(length(layout_widths)!=2){
        stop("Require 2 values in layout_widths when plotting all panels")
      }else{
        widths <- layout_widths
      }
    }
    if(is.null(layout_heights)){
      heights <- c(1,3)
    }else{
      if(length(layout_heights)!=2){
        stop("Require 2 values in layout_heights when plotting all panels")
      }else{
        heights <- layout_heights
      }
    }
    nf <- layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE), 
                 widths = widths, heights = heights)
    
  }else if(panels=="hist"){
    if(!is.null(layout_widths) || !is.null(layout_widths)){
      warning("layout_widths and layout_heights are ignored when plotting only 1 panel")
    }
    nf <- layout(matrix(1, ncol = 1), widths = 1, heights = 1)
    
  }else if(panels=="histforest"){
    if(!is.null(layout_widths)){
      warning("layout_widths are ignored when plotting only histogram + forest plot")
    }
    if(is.null(layout_heights)){
      heights <- c(1,3)
    }else{
      if(length(layout_heights)!=2){
        stop("Require 2 values in layout_heights when plotting all panels")
      }else{
        heights <- layout_heights
      }
    }
    nf <- layout(matrix(c(1, 2), ncol = 1, byrow = TRUE), 
                 widths = 1, heights = heights)
    
  }else if(panels=="forest"){
    if(!is.null(layout_widths) || !is.null(layout_widths)){
      warning("layout_widths and layout_heights are ignored when plotting only 1 panel")
    }
    nf <- layout(matrix(1, ncol = 1), widths = 1, heights = 1)
    
  }else if(panels=="foresttable"){
    if(is.null(layout_widths)){
      widths <- c(3,2)
    }else{
      if(length(layout_widths)!=2){
        stop("Require 2 values in layout_widths when plotting the forest plot with table")
      }else{
        widths <- layout_widths
      }
    }
    if(!is.null(layout_heights)){
      warning("layout_heights are ignored when plotting only the forest plot")
    }
    nf <- layout(matrix(c(1, 2), ncol = 2, byrow = TRUE), 
                 widths = widths, heights = 1)
    
  }else{
    stop("Bad argument to panels somehow missed by initial checks?")
  }
  
  
  
  ###
  # density histogram
  ###
  
  if(panels %in% c("all","hist","histforest")){
    
    par(mar=c(.1,3,.5,.1), mgp=c(1,.5,.0)) 
    
    hh <- hist(mc$ratio,
               breaks=c(-Inf,seq(min(forest$lb),max(forest$ub),length.out=num_breaks),Inf),
               plot = F)
    
    plot(hh$mids,hh$density,
         type="l",
         xlim=c(min(0, -0.1+min(forest$lb)),
                max(0, 0.1+max(forest$ub))),     
         ylab="Density",
         xlab="Effect size ratio",
         bty="l",
         cex.axis=0.75,
         cex.lab=0.75,
         xaxt='n',
         yaxt='n',
         yaxs='i')
    
    if(any(forest$cat==(k+1))){
      hh0 <- hist(mc$ratio[mc$id %in% forest$id[forest$cat==(k+1)]],
                  breaks=c(-Inf,seq(min(forest$lb),max(forest$ub),length.out=num_breaks),Inf),
                  plot = F)
      polygon(x=c(hh0$mids[2],hh0$mids[2:(length(hh0$density)-1)],tail(hh0$mids,2)[1]),
              y=mean(forest$cat==(k+1))*c(0,hh0$density[2:(length(hh0$density)-1)],0),
              col=adjustcolor(cols[k+1], alpha.f=0.7),
              border=NA)
    }
    
    for(i in 1:k){
      hhk <- hist(mc$ratio[mc$id %in% forest$id[forest$cat==i]],
                  breaks=c(-Inf,seq(min(forest$lb),max(forest$ub),length.out=num_breaks),Inf),
                  plot = F)
      polygon(x=c(hhk$mids[2],hhk$mids[2:(length(hhk$density)-1)],tail(hhk$mids,2)[1]),
              y=mean(forest$cat==i)*c(0,hhk$density[2:(length(hhk$density)-1)],0),
              col=adjustcolor(cols[i], alpha.f=0.7),
              border=NA)
    }
    
    
    for(i in 1:k){
      # separate loop to keep these reference lines on top layer of plot
      abline(v=selected_model$slopes[i,"b"],col=cols[i],lty=2,lwd=2)
    }
    
  }
  
  
  
  ###
  # scatter plot
  ###
  
  if(panels %in% c("all")){
    
    plot.magicmap(
      selected_model,
      which_model = NULL,
      class_thresh = class_thresh,
      diagnostics = NULL,
      colors = cols,
      legend = FALSE,
      conf_region = conf_region,
      hide_warnings = TRUE, # already warns outside this outside function
      mar=c(1,1,.5,.1),
      mgp=c(0,0,0),
      xlim=c(0,1.1*max(tab[,1]+1.75*tab[,2])),
      ylim=c(min(0,min(tab[,3]-1.75*tab[,4])),
             max(0,max(tab[,3]+1.75*tab[,4]))),
      cex.axis=0.3, 
      cex.lab=0.6,
      xlab=paste(target_name,"GWAS Beta"), 
      ylab=paste(comparator_name,"GWAS Beta"), 
      xaxt='n', 
      yaxt='n',
      bty='n',
      axes=F,
      xaxs="i",
      # yaxs="i",
      lwd=2,
      se_length=0.02,
      se_color="gray85",
      assigned_cex=1.2, 
      unassigned_cex=0.75
    )
    
  }
  
  
  
  ###
  # forest plot
  ###
  
  if(panels %in% c("all","histforest","forest","foresttable")){
    
    par(mar=c(2,3,.1,.1), mgp=c(1,0.5,0)) 
    
    # sort by precision and optionally drop unassigned obs
    idx_tmp <- order(forest$ub-forest$lb,decreasing = T)
    idx <- NULL
    if(forest_drop_unassigned){
      forest_comps <- c(k:1)
    }else{
      forest_comps <- c((k+1):1)
    }
    for(cc in forest_comps){
      idx <- c(idx, idx_tmp[forest$cat[idx_tmp]==cc])
    }
    forest_include <- forest[idx,]
    
    # init plot
    plot(c(0,2),
         col=rgb(0,0,0,0),
         xlim=c(min(0, -0.1+min(forest$lb)),
                max(0, 0.1+max(forest$ub))),
         ylim=c(1,nrow(forest_include)+2),
         bty='n',
         xlab="Effect size ratio",
         ylab="",
         cex.axis=0.5,
         cex.lab=0.5,
         yaxt='n'
    )
    
    
    for(i in 1:nrow(forest_include)){
      arrows(x0=forest_include$lb[i],x1=forest_include$ub[i],
             y0=i,y1=i,
             code=3,
             angle = 90,
             col = "grey80",
             length=0.02
      )
    }
    
    points(forest_include$rat,1:nrow(forest_include),
           pch=20,
           cex=ifelse(forest_include$cat==0, 1, 1.5),
           col=forest_include$col)
    
    mtext(substitute(paste(bold("Locus"))),
          side=2,
          at=nrow(forest_include)+1.25,
          las=2,
          cex=1.33*cex_forest_labels,
          adj=1.25)
    
    for(i in 1:k){
      abline(v=selected_model$slopes[i,"b"],col=cols[i],lty=2,lwd=2)
    }
    
    axis(2, 1:nrow(forest_include), labels = forest_include$id, tick = F, las=2, cex.axis=cex_forest_labels)
    
  }
  
  
  
  ###
  # forest plot table
  ###
  
  if(panels %in% c("all","foresttable")){
    
    par(mar=c(2,.1,.1,.1), mgp=c(1.5,0.5,0)) # c(5,4,4,2)+.1
    
    plot(c(0,2),
         col=rgb(0,0,0,0),
         xlim=c(-1,1),
         ylim=c(1,nrow(forest_include)+2),
         bty='n',
         xlab="",
         ylab="",
         xaxt='n',
         cex.axis=0.6,
         cex.lab=0.75,
         yaxt='n'
    )
    
    text(paste0(sprintf("%0.3f",tab[idx,1])," (",sprintf("%0.3f",tab[idx,2]),")"),x=-.8,y=1:nrow(forest_include),cex=cex_forest_labels)
    text(paste0(sprintf("%0.3f",tab[idx,3])," (",sprintf("%0.3f",tab[idx,4]),")"),x=0,y=1:nrow(forest_include),cex=cex_forest_labels)
    text(sprintf("%0.3f",forest_include$prob),x=.8,y=1:nrow(forest_include),cex=cex_forest_labels)
    text(x=-.8,y=nrow(forest_include)+2,adj=0.5,bquote(paste(bold(.(target_name)))), cex=1.33*cex_forest_labels)
    text(x=-.8,y=nrow(forest_include)+1,adj=0.5,bquote(paste(bold("Beta (SE)"))), cex=1.33*cex_forest_labels)
    text(x=0,y=nrow(forest_include)+2,adj=0.5,bquote(paste(bold(.(comparator_name)))), cex=1.33*cex_forest_labels)
    text(x=0,y=nrow(forest_include)+1,adj=0.5,bquote(paste(bold("Beta (SE)"))), cex=1.33*cex_forest_labels)
    text(x=.8,y=nrow(forest_include)+2,adj=0.5,bquote(paste(bold("Posterior"))), cex=1.33*cex_forest_labels)
    text(x=.8,y=nrow(forest_include)+1,adj=0.5,bquote(paste(bold("Prob."))), cex=1.33*cex_forest_labels)
    
  }
  
  nf <- layout(matrix(1, ncol = 1), widths = 1, heights = 1)
  par(par_backup)
  
  
  ###
  # optional return values to save data
  ###
  if(output_forest_data=="raw"){
    return(mc)
  }else if(output_forest_data=="forest"){
    return(forest)
  }
  
}


