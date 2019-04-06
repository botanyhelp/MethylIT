#' @rdname print.glmDataSet
#' @aliases print.glmDataSet
#' @aliases print.RangedGlmDataSet
#' @aliases print
#' @title Printing object from \emph{glmDataSet} and \emph{"RangedGlmDataSet"}
#'     classes by simple print methods.
#' @param x Object from class \emph{glmDataSet} or from class
#'     \emph{"RangedGlmDataSet"}.
#' @param digits Number of significant digits to be used.


#' @rdname print.glmDataSet
#' @name print.glmDataSet
#' @export
#' @keywords internal
print.glmDataSet <- function(x, ...) UseMethod("print", x)
print.glmDataSet <- function(x, digits = getOption("digits")) {
   cm <- dim(counts)
   lvs <- levels(x$colData$condition)
   cat("glmDataSet with ", cm[1], " rows and ", cm[2], " columns (individuals)",
       " with factor levels ", "'", lvs[1],"'", " and '", lvs[2],"' \n",
       sep = "")
   cat("The accessible objects in the dataset are: \n")
   print(summary(x), digits = digits)
   cat("\n")
}

#' @rdname print.glmDataSet
#' @name print.RangedGlmDataSet
#' @export
#' @keywords internal
print.RangedGlmDataSet <- function(x, ...) UseMethod("print", x)
print.RangedGlmDataSet <- function(x, digits = getOption("digits")) {
    r <- length(x$GR)
    lvs <- levels(x$colData$condition)
    col <- length(x$colData$condition)
    cat("RangedGlmDataSet with ", r, " regions and ", col,
       " columns (individuals) ",
       "with factor levels ", "'", lvs[1],"'", " and '", lvs[2],"' \n",
       sep = "")
  cat("The accessible objects in the dataset are: \n")
  print(summary(x))
  cat("\n")
}

