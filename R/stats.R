


##'
calculatePCA <-function(){
  
}

##' 
performQC <-function(){
  
}


##' QQ Plot for P-values
##'
##' This function generates a Quantile-Quantile (QQ) plot for a given set of p-values. It also allows 
##' for optional thresholds based on Bonferroni correction and False Discovery Rate (FDR) adjustment, 
##' as well as a user-defined line.
##'
##' @param Pval A numeric vector of p-values.
##' @param bonf Logical; if TRUE, applies a Bonferroni correction and adds a threshold line to the plot. Default is TRUE.
##' @param fdr.bh Logical; if TRUE, applies an FDR (Benjamini-Hochberg) correction and adds a threshold line to the plot. Default is TRUE.
##' @param userwideline A numeric value to specify a custom y-intercept for an additional horizontal line on the plot. Default is NULL.
##' @param thresh A numeric value specifying the significance threshold for the Bonferroni and FDR corrections. Default is 0.05.
##'
##' @return A ggplot2 object representing the QQ plot with optional threshold lines.
##' @author Jemay SALOMON
##' @examples
##' # Example usage:
##' p_values <- runif(1000, min = 0, max = 1)
##' qqPlot(p_values)
##'
##' @export
qqplotGWAS <- function(Pval, bonf = TRUE, fdr.bh = TRUE,
                       userwideline = NULL, thresh = 0.05,
                       verbose = FALSE, plot.nb.pval = FALSE) {

  packageCheck("ggplot2")
  Pval <- Pval[!is.na(Pval)]
  stopifnot(is.numeric(Pval), is.logical(bonf),
            is.logical(fdr.bh))
  if(!is.null(userwideline))
    stopifnot(is.numeric(userwideline))
  
  
  #check NA
  #add confidence interval
  
  observed <- -log10(sort(Pval))
  expected <- -log10(ppoints(length(Pval)))
  
  bonf <- if (bonf) -log10(thresh / length(Pval)) else NULL
  fdr.bh <- if (fdr.bh) {
    pv_bh <- stats::p.adjust(Pval, method = "BH")
    -log10(sort(Pval[pv_bh <= thresh], decreasing = TRUE)[1])
  } else NULL
  
  xdat <- data.frame(Expected = expected, Observed = observed)
  maxlogP <- max(c(observed, expected, bonf, fdr.bh, userwideline), na.rm = TRUE)
  
  p <- ggplot(xdat, aes(x = Expected, y = Observed)) +
    geom_point(size = 1.3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = "QQ Plot") +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA), 
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),     
      panel.grid.major = element_blank(),                     
      panel.grid.minor = element_blank()  
    ) +
    xlim(0, maxlogP) +
    ylim(0, maxlogP)
  
  if (!is.null(bonf)) {
    p <- p + geom_hline(yintercept = bonf, color = "blue",
                        linetype = "solid", linewidth = 0.5) +
      annotate("text", x = max(maxlogP), y = bonf, 
               label = "bonf", vjust = -0.5, color = "blue")
  }
  
  if (!is.null(fdr.bh)) {
    p <- p + geom_hline(yintercept = fdr.bh, color = "green", 
                        linetype = "solid", linewidth = 0.5) +
      annotate("text", x = max(maxlogP), y = fdr.bh, 
               label = "fdr.bh", vjust = -0.5, color = "green")
  }
  
  if (!is.null(userwideline)) {
    p <- p + geom_hline(yintercept = userwideline, color = "black", 
                        linetype = "solid", linewidth = 0.5) +
      annotate("text", x = max(maxlogP), y = userwideline, 
               label = "userwideline", vjust = -.5, color = "black")
  }
  
  print(p)
}


##'
annotateVariants <-function(){
  
}


##'
plotPCA <-function(){
  
}


##'
summaryReport <-function(){
  
}

##' Calculates polygenic risk scores (PRS) using GWAS results and test data.
polygenicRiskScore <- function(){
  
}


##' Impute Missing Genotypes
##'
##' This function imputes missing genotype data in a numeric matrix using either the mean or mode. The imputation is performed row-wise, with genotypes in rows and SNPs in columns.
##'
##' @param X A numeric matrix with genotypes in rows and SNPs in columns, with missing values represented as `NA`.
##' @param with A character string specifying the imputation method. Either `"mean"` or `"mode"`. Default is `"mean"`.
##' @return A numeric matrix with missing values imputed.
##' @author Jemay SALOMON
##' @examples
##' # Example usage:
##' X <- matrix(c(0, 1, NA, 2, 2, NA, 1, 0, 1, NA), nrow = 5, byrow = TRUE)
##' rownames(X) <- c("Ind1", "Ind2", "Ind3", "Ind4", "Ind5")
##' colnames(X) <- c("SNP1", "SNP2")
##' imputedX <- imputeMissingGenotypes(X, with = "mean")
##' print(imputedX)
##' imputedX <- imputeMissingGenotypes(X, with = "mode")
##' print(imputedX)
##' 
##' @export
imputeMissingGenotypes <- function(X, with = "mean") {
  
  # Ensure X is a numeric matrix
  stopifnot(is.matrix(X), is.numeric(X))
  
  # Check if method is valid
  if (!with %in% c("mean", "mode")) {
    stop("Method must be either 'mean' or 'mode'.")
  }
  
  # Check if there are any missing values
  if (!any(is.na(X))) {
    stop("No NAs to impute.")
  }
  
  if(with=="mode"){
    X <- apply(X, 2, function(x){
      freq <- table(x)
      x[is.na(x)] <- as.integer(names(which.max(freq)))
      return(x)
    })
  }
  
  if(with=="mean"){
    
    #TODO :implement
    
  }
  
  return(X)
}


##' 
geneSetEnrichment<-function(){
}


##' 
filterVariants<-function(){
}