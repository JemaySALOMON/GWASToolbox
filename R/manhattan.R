
##' Manhattan Plot Function
##'
##' Creates a Manhattan plot for genome-wide association studies (GWAS) data. 
##' The plot visualizes SNPs across chromosomes, highlighting those that are 
##' significant based on Bonferroni correction, FDR adjustment, and optionally 
##' a user-defined threshold. SNPs of interest can be annotated with labels.
##'
##' @param X A data frame containing the GWAS data.
##' @param CHR The column name (as a string) in `X` that contains chromosome numbers.
##' @param BP The column name (as a string) in `X` that contains base pair positions.
##' @param Pval The column name (as a string) in `X` that contains p-values.
##' @param SNP The column name (as a string) in `X` that contains SNP identifiers.
##' @param bonf Logical value indicating if Bonferroni correction should be applied.
##' @param fdr.bh Logical value indicating if FDR adjustment should be applied.
##' @param thresh The threshold for significance testing.
##' @param userwideline Optional numeric value specifying an additional horizontal line 
##'   to be added to the plot. If not NULL, it will be used to draw a line on the plot.
##' @param show.labels Optional vector of SNPs to be highlighted and labeled on the plot.
##'
##' @return A ggplot object of the Manhattan plot.
##' @import ggplot2 dplyr ggrepel
##'
##' @author Jemay Salomon
##' @examples
##' # Example usage
##' manhattanPlot(X = my_data, CHR = "chr", BP = "bp", Pval = "pval", SNP = "snp", 
##'                bonf = TRUE, fdr.bh = TRUE, thresh = 0.05, 
##'                userwideline = NULL, show.labels = c("snp001", "snp002"))
##' 
##' @export
manhattanPlot<-function(X, CHR = NULL, BP, Pval, SNP, bonf = TRUE, 
                        fdr.bh = TRUE,thresh = 0.05, userwideline = NULL,
                        show.labels = NULL, chrLabs = NULL) {
  
  # Load necessary packages
  packageCheck("ggrepel")
  packageCheck("ggplot2")
  packageCheck("dplyr")
  
  
  # Check inputs for existence and type
  stopifnot(
    BP %in% names(X),
    Pval %in% names(X),
    SNP %in% names(X),
    is.numeric(X[[BP]]),
    is.numeric(X[[Pval]]),
    is.numeric(thresh),
    is.logical(bonf),
    is.logical(fdr.bh)
  )
  
  # Extract dynamic cols
  SNP <-X[[SNP]]
  P <- X[[Pval]]
  
  
  if (!is.null(CHR)) {
    stopifnot(CHR %in% names(X),
              is.numeric(X[[CHR]]))
  }
  
  if (!is.null(chrLabs)) {
    stopifnot(chrLabs %in% names(X))
  }
  
  if (!is.null(userwideline)) {
    stopifnot(is.numeric(userwideline))
  }
  
  # Ensure either CHR or chrLabs is provided
  if (is.null(CHR) && is.null(chrLabs)) {
    stop("Either 'CHR' or 'chrLabs' must be provided and present in the data frame.")
  }
  
  # Determine which column to use for chromosomes
  chromosome_col <- if (!is.null(chrLabs)) chrLabs else CHR
  if (!chromosome_col %in% names(X)) {
    stop(paste("The column", chromosome_col, "is not found in the data frame."))
  }
  
  # Threshold calculations
  bonf <- if (bonf) thresh / length(X[[Pval]]) else NULL
  
  fdr.bh <- if (fdr.bh) {
    pv.bh <- stats::p.adjust(X[[Pval]], method = "BH")
    sort(X[[Pval]][pv.bh <= thresh], decreasing = TRUE)[1]
  } else NULL
  
  # Convert the column name to a symbol for dplyr functions
  chromosome_col_sym <- rlang::sym(chromosome_col)
  
  x <- X %>%
    group_by(!!chromosome_col_sym) %>%
    summarise(chr_len = max(!!sym(BP)), .groups = 'drop') %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(X, ., by = setNames(chromosome_col, chromosome_col)) %>%
    arrange(!!chromosome_col_sym, !!sym(BP)) %>%
    mutate(cumBP = !!sym(BP) + tot)
  
  
  # Highlight SNPS of interest
  if (!is.null(show.labels)) {
    if (!is.null(bonf) && !is.null(fdr.bh)) {
      x <- x %>%
        mutate(
          is_highlight = ifelse(SNP %in% show.labels, "yes", "no"),
          is_annotate = ifelse(
            is_highlight == "yes" & (
              -log10(P) > -log10(bonf) | -log10(P) > -log10(fdr.bh)
            ), 
            "yes", 
            "no"
          )
        )
    } else if (!is.null(bonf)) {
      x <- x %>%
        mutate(
          is_highlight = ifelse(SNP %in% show.labels, "yes", "no"),
          is_annotate = ifelse(
            is_highlight == "yes" & (
              -log10(P) > -log10(bonf)
            ), 
            "yes", 
            "no"
          )
        )
    } else if (!is.null(fdr.bh)) {
      x <- x %>%
        mutate(
          is_highlight = ifelse(SNP %in% show.labels, "yes", "no"),
          is_annotate = ifelse(
            is_highlight == "yes" & (
              -log10(P) > -log10(fdr.bh)
            ), 
            "yes", 
            "no"
          )
        )
    } else {
      x <- x %>%
        mutate(
          is_highlight = ifelse(SNP %in% show.labels, "yes", "no"),
          is_annotate = ifelse(is_highlight == "yes", "yes", "no")
        )
    }
  }
  
  # Highlight significant SNPs
  if (is.null(show.labels)) {
    if (!is.null(fdr.bh) &&  is.null(bonf)) {
      x <- x %>%
        mutate(
          is_highlight = ifelse(-log10(P) > -log10(fdr.bh), "yes", "no")
        )
    } else if (is.null(fdr.bh) && !is.null(bonf)) {
      x <- x %>%
        mutate(
          is_highlight = ifelse(-log10(P) > -log10(bonf), "yes", "no")
        )
    }else if (!is.null(fdr.bh) && !is.null(bonf)) {
      x <- x %>%
        mutate(
          is_highlight = ifelse(
            -log10(P) > -log10(bonf) | -log10(P) > -log10(fdr.bh), 
            "yes", 
            "no"
          )
        )
    } else {
      x <- x %>%
        mutate(
          is_highlight = "no"
        )
    }
  }
  
  # Prepare axis data
  axisdat <- x %>%
    group_by(!!chromosome_col_sym) %>%
    summarize(idxcenter = (max(cumBP) + min(cumBP)) / 2)
  
  # Point colors
  point.colors <- rep(c("grey", "skyblue"), length.out = length(unique(x[[chromosome_col]])))
  names(point.colors) <- unique(x[[chromosome_col]])
  
  # Plot
  p <- ggplot(x, aes(x = cumBP, y = -log10(P))) +
    geom_point(aes(color = as.character(!!chromosome_col_sym)), alpha = 0.8, size = 1.3) +
    scale_color_manual(
      values = c(
        "bonf" = "blue",
        "fdr.bh" = "green",
        "userwideline" = "black",
        point.colors
      ),
      name = "Correction",
      breaks = c("bonf", "fdr.bh", "userwideline")
    ) +
    theme_bw() +
    scale_x_continuous(labels = axisdat[[chromosome_col_sym]], breaks = axisdat$idxcenter) +
    theme(
      legend.position = "right",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Chromosome", y = "-log10(P)")
  
  # Add horizontal lines
  if (!is.null(userwideline))
    p <- p + geom_hline(aes(yintercept = userwideline, color = "userwideline"), linetype = "solid", linewidth = 0.5)
  if (!is.null(bonf))
    p <- p + geom_hline(aes(yintercept = -log10(bonf), color = "bonf"), linetype = "solid", linewidth = 0.5)
  if (!is.null(fdr.bh))
    p <- p + geom_hline(aes(yintercept = -log10(fdr.bh), color = "fdr.bh"), linetype = "solid", linewidth = 0.5)
  
  # Add highlighted SNPs and labels
  if (!is.null(show.labels)) {
    p <- p +
      geom_point(data = subset(x, is_highlight == "yes" & is_annotate == "yes"), color = "orange", size = 2) +
      geom_label_repel(
        data = subset(x, is_highlight == "yes" & is_annotate == "yes"),
        aes(label = SNP),
        size = 2,
        max.overlaps = Inf
      )
  } else {
    p <- p +
      geom_point(data = subset(x, is_highlight == "yes"), color = "orange", size = 2)
  }
  
  print(p)
}


manhattanPlot(dat, CHR = "chr", BP ="bp", Pval="P", SNP="snp", bonf = TRUE, 
                        fdr.bh = T,thresh = 0.05, userwideline = 6,
                        show.labels = NULL, chrLabs = NULL)
