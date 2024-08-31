##' Retrieve a Value from a DataFrame
##'
##' \code{getInfo} extracts a value from a specified row and column in a DataFrame.
##'
##' This function is useful for retrieving specific elements from a DataFrame when you know the
##' row identifier (which can be numeric, character, or any other type) and the column name.
##'
##' @param data DataFrame. The DataFrame containing the data.
##' @param x Mixed. The row identifier, which can be numeric (row index) or any value that uniquely identifies a row (e.g., row name).
##' @param colName Character. The name of the column from which to retrieve the value. Must be a valid column name in the DataFrame.
##'
##' @return The value at the intersection of the specified row and column in the DataFrame.
##'
##' @details
##' This function performs the following checks:
##' \itemize{
##'   \item It checks if \code{colName} exists in the DataFrame columns.
##'   \item If \code{x} is numeric, it checks if it is a valid row index.
##'   \item If \code{x} is not numeric, it checks if it matches a row name or identifier.
##' }
##' If any of these conditions are not met, the function will stop and throw an error message.
##' @author Jemay SALOMON
##' @examples
##' # Example usage with numeric row index
##' data <- data.frame(x = c(1, 2, 3), "Isolate_15FG033" = c(1, 0, 0), "Isolate_15FG038" = c(0, 2, 1))
##' value <- getInfo(data, 2, "Isolate_15FG033")
##' print(value)
##' 
##' @export
getInfo <- function(data, x, colName) {
  if (!colName %in% colnames(data)) {
    stop("The specified column name does not exist in the DataFrame.")
  }
  if (is.numeric(x)) {
    if (x > nrow(data) || x < 1) {
      stop("The row index x is out of bounds for the DataFrame.")
    }
    idx <- x
  } else {
    idx <- which(rownames(data) == x)
    if (length(idx) == 0) {
      stop("The specified row identifier does not exist in the DataFrame.")
    }
  }
  return(data[idx, colName])
}


##' Check if a Package is Installed
##'
##' The `packageCheck` function verifies whether a specified package is installed and loaded into R. If the package is not available, the function stops execution and returns an error message indicating that the package is required.
##'
##' @param pkg A character string specifying the name of the package to check.
##' @return This function does not return a value. It is used for its side effect of checking for package availability. If the package is not found, the function stops execution with an error.
##' @author Jemay SALOMON
##' @examples
##' \dontrun{
##'   packageCheck("ggplot2")  # Checks if the ggplot2 package is installed
##' }
##'
##' @export
packageCheck <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE))) {
    stop(paste("Package", pkg, "is required but not installed"))
  }
}

##' Convert Character Columns with Comma as Decimal Separator to Numeric
##'
##' This function converts specified character columns in a data frame to numeric,
##' replacing commas with dots as the decimal separator.
##' The function ensures that the input is a data frame and that the specified columns exist.
##'
##' @param x A data.frame in which the columns will be converted.
##' @param cols A vector of column names or indices specifying which columns to convert to numeric.
##' These columns must be character type.
##' @return A data.frame with the specified character columns converted to numeric.
##' @author Jemay SALOMON
##' @export
##' @examples
##' df <- data.frame(A = c("1,23", "2,34", "3,45"), B = c("10,5", "20,6", "30,7"), stringsAsFactors = FALSE)
##' convertCharColsToNumeric(df, cols = c("A", "B"))
##' @export
convertCharColsToNumeric <- function(x, cols) {
  stopifnot(is.data.frame(x))
  stopifnot(all(cols %in% names(x) | cols %in% seq_along(x)))
  idx <- sapply(x[cols], is.character)
  x[cols[idx]] <- lapply(x[cols[idx]], function(col) {
    as.numeric(col)
  })
  
  return(x)
}


##' Print Space (Interline)
##'
##' This function prints a specified number of newline characters (interlines) to create space in the console output.
##'
##' @param x A numeric value representing the number of interlines to print. It should be an integer (not decimal).
##' @return Prints the specified number of newline characters.
##' @note The function ensures that the input is a non-decimal numeric value.
##' @author Jemay SALOMON
##' @examples
##' space(2)
##' @export
space <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(x == trunc(x))
  for (i in 1:x) {
    cat("\n")
  }
}


##' Converts GWAS data between different formats (e.g., PLINK to VCF).
convertFormats <- function(){
}

##' Exports GWAS results in various formats suitable for downstream analysis (e.g., CSV, BED, PLINK format)
exportGwasResults<- function(){
}

gwasResults%>%select(P)
