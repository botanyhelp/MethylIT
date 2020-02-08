#' @rdname glmDataSet
#' @name glmDataSet
#' @title Data set constructor for class glmDataSet
#' @description This function is used to build a object suitable to be used
#'     with Methyl-IT \code{link{countTest2}} function.
#' @details Data set constructor for class glmDataSet also validate the object
#' @param GR A GRanges object with the count matrix of DMPs in the metacolumns
#'     (see \emph{'counts'}). If provided, then leave parameter
#'     \emph{'counts = NULL'}.
#' @param counts Count matrix of DMPs with minimal dimensions 1 (row) x 4
#'     (columns). Column names must corresponds to the rownames from parameter
#'      'colData'.
#' @param colData A data frame with one columnn named 'condition', which must be
#'     a factor with exactly two levels. The rownames of \emph{colData}
#'     individual samples. The row names of \emph{colData} must correspond to
#'     th column names of the count matrix.
#' @export
#' @examples
#' set.seed(133) # Set a seed
#' ## A GRanges object with the count matrix in the metacolumns is created
#' countData <- matrix(sample.int(200, 500, replace = TRUE), ncol = 4)
#' colnames(countData) <- c("A1","A2","B1","B2")
#' start <- seq(1, 25e4, 2000)
#' end <- start + 1000
#' chr <- c(rep("chr1", 70), rep("chr2", 55))
#' GR <- GRanges(seqnames = chr, IRanges(start = start, end = end))
#' mcols(GR) <- countData
#' ## Gene IDs
#' names(GR) <- paste0("gene", 1:length(GR))
#'
#' ## An experiment design is set.
#' colData <- data.frame(condition = factor(c("A","A","B","B")),
#'                       c("A1","A2","B1","B2"),
#'                       row.names = 2)
#' ## A RangedGlmDataSet is created
#' ds <- glmDataSet(GR = GR, colData = colData)
#'
glmDataSet <- function(GR = NULL, counts = NULL, colData = NULL) {
   if (is.null(GR) && is.null(counts)) {
       cat("\n")
       stop("'GR' or 'counts' must be provided")
   }
   if (is.null(colData$condition)) {
      cat("\n")
      stop("In 'colData', 'condition' must be provided")
   }
    if (!is(colData$condition, "factor")) {
       cat("\n")
       stop("In 'colData','condition' must be a 'factor'")
   }
   if (!is.null(colData$condition)) {
      if (length(levels(colData$condition)) != 2) {
         cat("\n")
         stop("In 'colData', factor 'condition' must have only two levels")
      }
   }
   if (!is.null(GR) && !is(GR, "GRanges")) {
       cat("\n")
       stop("'GR' must be a GRanges")
   }
   if (!is.null(counts) && !is(counts, "matrix")) {
       counts <- try(as.matrix(counts), silent = TRUE)
       if (inherits(counts, "try-error")) {
           cat("\n")
           stop("'counts' must be a matrix or an object coercible as a matrix")
       }
   }
   if (!is.null(counts)) {
       if (is.null(colnames(counts))) {
           cat("\n")
           stop("The 'counts' matrix must have column names, which \n",
                "must correspond to 'colData' rownames")
       }
       if (any(colnames(counts) != rownames(colData))) {
           cat("\n")
           stop("Individual names in 'counts' must correspond to ",
              "'colData' rownames")
       }
   }
   if (is.null(counts)) {
       if (ncol(mcols(GR)) > 0) {
           if (any(rownames(colData) != colnames(mcols(GR)))) {
               cat("\n")
               stop("Individual names in 'GR' must correspond to \n",
                   "'colData' rownames")
           }
       } else stop("Since 'counts = NULL', GR must have a count matrix in the",
                   " metacolumns or provide 'counts'")
   }
   if (is.null(GR)) {
       if (any(rownames(colData) != colnames(counts))) {
           cat("\n")
           stop("Individual names in 'counts' must correspond to \n",
               "condition rownames")
       }
   }

   if (is.null(GR)) {
       x <- list(counts = counts, colData = colData,
                 sampleNames = rownames(colData),
                 levels = levels(colData$condition),
                 optionData = NULL)
       x <- structure(x, class = c("glmDataSet"))
   }
   else {
       if (is.null(counts)) {
           counts <- as.matrix(mcols(GR))
           mcols(GR) <- NULL
       }
       # To make sure information is not redundant
       if (!is.null(counts) && !is.null(GR)) mcols(GR) <- NULL

       # The final output
       x <- list(GR = GR, counts = counts, colData = colData,
                 sampleNames = rownames(colData),
                 levels = levels(colData$condition),
                 optionData = NULL)
       x <- structure(x, class = c("RangedGlmDataSet"))
   }
   return(x)
}


