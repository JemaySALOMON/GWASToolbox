##' Run GWAS using MM4LMM, mlmm, eagle Methods
##'
##' The `gwasRun` function performs a Genome-Wide Association Study (GWAS) using one of three specified methods: MM4LMM, mlmm, or eagle. The function takes phenotype data, genotype data, and variance-covariance matrices as input and returns a list containing the results of the analysis.
##'
##' @param Y A numeric vector representing the phenotype data.
##' @param listXs A list containing the genotype data matrix, with the key `"X"` representing the matrix of genotypes.
##' @param listVCov A list containing variance-covariance matrices. The key `"K"` represents the kinship matrix, and optionally, the key `"Error"` can represent an error variance matrix.
##' @param maxsteps An optional numeric value specifying the maximum number of steps for the `mlmm` method. This argument is required if the `method` is `"mlmm"`.
##' @param method A character string specifying the method to use for the GWAS. Options are `"mm4lmm"`, `"mlmm"`, or `"eagle"`.
##' @param verbose A logical value indicating whether to print progress messages. Default is `TRUE`.
##'
##' @return A list containing the results of the GWAS. The structure of the list varies depending on the method used:
##' - For `"mm4lmm"`: The list includes the estimated model (`mmest`), test results (`out.test`), and p-values (`res.mm4lmm`).
##' - For `"mlmm"`: The list includes the results of the `mlmm` function (`mygwas`).
##' - For `"eagle"`: The list includes the results of the `eagle`
##'
##' @details 
##' - **MM4LMM**: Utilizes the `MM4LMM` package to fit a linear mixed model and compute p-values.
##' - **mlmm**: Uses the `mlmm` package to perform a multi-locus mixed-model analysis.
##' - **eagle**: Intended for use with the `eagle` package to perform GWAS
##'
##' @author Jemay SALOMON
##' @examples
##' \dontrun{
##'   # Example usage for the mm4lmm method
##'   result <- gwasRun(Y = phenotype_data, listXs = genotype_list, listVCov = covariance_list, method = "mm4lmm")
##'
##'   # Example usage for the mlmm method
##'   result <- gwasRun(Y = phenotype_data, listXs = genotype_list, listVCov = covariance_list, maxsteps = 6, method = "mlmm")
##'   
##'   # Example usage for the eagle method
##'   
##' }
##'
##' @export
gwasRun <- function(Y, listXs, listVCov,  
                    maxsteps = 6, method=NULL, 
                    verbose=TRUE){
  
  # range-check
  stopifnot(is.numeric(Y), 
            !is.null(listXs),
            !is.null(listVCov),
            is.numeric(listXs[["X"]]))
  stopifnot(length(Y) == nrow(listVCov[["K"]]))
  stopifnot(nrow(listVCov["K"]) == nrow(listXs[["X"]]))
  
  # conditions
  switch(method,
         
         "mm4lmm" = { 
           if(verbose)
             print("check packages MM4LMM")
           packageCheck("MM4LMM")
           
           if(verbose)
             print("fit MM4LMM")
           Error = if (is.null(listVCov$Error)) diag(length(Y)) else listVCov$Error
           VarList = list(Additive = listVCov[["K"]] , Error = Error); rm(Error)
           mmest <- MMEst(Y = Y, X = listXs[["X"]], VarList = VarList)
           
           if(verbose)
             print("compute p-value [MM4LMM]")
           out <- list(
             mmest = mmest,
             out.test = AnovaTest(mmest),
             res.mm4lmm = cbind(P = sapply(AnovaTest(mmest), function(x) if("Xeffect" %in% rownames(x)) {
               return(x["Xeffect", "pval"])
             } else {
               return(NA)
             }))
           )
         },
         
         "mlmm" = {
           if(verbose)
             print("check packages mlmm")
           packageCheck("mlmm")
           stopifnot(!is.null(maxsteps))
           
           if(verbose)
             print("fit mlmm")
           out <- list(
             mlmm = mlmm::mlmm(Y = Y, X = listXs[["X"]], K = listVCov[["K"]], maxsteps = maxsteps, nbchunks = 2)
           )
         },
         
         "eagle" = {
           if(verbose)
             print("check packages eagle")
           packageCheck("eagle")
           stopifnot(!is.null(maxsteps))
           
           if(verbose)
             print("fit eagle")
           out <- list(
            #TODO : implement eagle
           )
         },
         
         stop("method must be one of these: c('mm4lmm','mlmm','eagle')")
  )
  return(out)
}

##' VCModel Class
##'
##' The `VCModel` class is designed to extract variance components (`VA`, `VE`, and `h2`) from different genetic analysis methods.
##' Currently, it supports the `mm4lmm` method, with placeholders for `mlmm` and `eagle` methods.
##'
##' @field method A character string indicating the method used for variance component extraction (e.g., "mm4lmm").
##' @field obj A list containing the results from the variance component estimation method.
##' @field VA A numeric value representing the additive genetic variance component.
##' @field VE A numeric value representing the residual (environmental) variance component.
##' @field h2 A numeric value representing the heritability estimate (VA / (VA + VE)).
##' @author Jemay SALOMON
##' @export VCModel
##' @examples
##' # Example for using VCModel with the MM4LMM method
##' mmest <- list(list(Sigma2 = c(2.0, 3.0)))
##' vc <- VCModel$new(obj = mmest, method = "mm4lmm")
##' vc$show()

VCModel <- setRefClass(
  "VCModel",
  fields = list(
    method = "character",
    obj = "list",
    VA = "numeric",
    VE = "numeric",
    h2 = "numeric"
  ),
  
  methods = list(
    initialize = function(obj, method = "mm4lmm") {
      stopifnot(
        is.list(obj),
        length(obj) > 0,
        is.character(method)
      )
      
      .self$obj <- obj
      .self$method <- method
      
      if (method == "mm4lmm") {
        .self$initializeMM4LMM()
      } else if (method == "mlmm") {
        .self$initializeMLMM()
      } else if (method == "eagle") {
        .self$initializeEagle()
      } else {
        stop("Unknown method specified.")
      }
    },
    
    initializeMM4LMM = function() {
      stopifnot(
        is.list(.self$obj[[1]]),
        "Sigma2" %in% names(.self$obj[[1]]),
        length(.self$obj[[1]]$Sigma2) == 2
      )
      
      .self$VA <- .self$obj[[1]]$Sigma2[1]
      .self$VE <- .self$obj[[1]]$Sigma2[2]
      .self$h2 <- .self$VA / (.self$VA + .self$VE)
    },
    
    initializeMLMM = function() {
      stop("MLMM method is not yet implemented.")
    },
    
    initializeEagle = function() {
      stop("Eagle method is not yet implemented.")
    },
    
    getVA = function() {
      return(.self$VA)
    },
    
    getVE = function() {
      return(.self$VE)
    },
    
    getH2 = function() {
      return(.self$h2)
    },
    
    show = function() {
      cat("Method:", .self$method, "\n")
      cat("VA:", .self$VA, "\n")
      cat("VE:", .self$VE, "\n")
      cat("h2:", .self$h2, "\n")
    }
  )
)
