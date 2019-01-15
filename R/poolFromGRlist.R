#' @rdname poolFromGRlist
#'
#' @title Methylation pool from a list of GRanges objects with methylation read
#'     counts
#' @description This function will build a GRanges methylation pool from a list
#'     of GRanges objects
#' @details The list of GRanges objects (LR) provided to build a virtual
#'     methylome should be an output of the function 'readCounts2GRangesList' or
#'     at least each GRanges must have the columns named "mC" and "uC", for the
#'     read counts of methylated and unmethylated cytosines, respectively.
#'
#' @param LR  list of GRanges objects to build a virtual individual (methylation
#'     pool)
#' @param stat statistic used to estimate the methylation pool: row sum, row
#'     mean or row median of methylated and unmethylated read counts across
#'     individuals
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see bplapply function from
#'     BiocParallel package).
#' @param tasks integer(1). The number of tasks per job. Value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#' @param prob Logic. Whether the variable for pooling is between 0 and 1 (a
#'     probability), e.g., methylation levels. If TRUE, then Fisher's
#'     transformation is applied, the row mean is computed for each cytosine
#'     site and returned in the original measurement scale between 0 and 1 by
#'     using the inverse of Fisher's transformation.
#' @param column If prob == TRUE, then the 'column' from the LR metacolumns
#'     where the prob values are found must be provided. Otherwise, column = 1L.
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... Additional parameters for 'uniqueGRanges' function.
#'
#' @return A GRanges object
#'
#' @examples
#' gr1 <- makeGRangesFromDataFrame(
#'     data.frame(chr = "chr1", start = 11:15, end = 11:15,
#'                strand = '*', mC = 1, uC = 1:5),
#'     keep.extra.columns = TRUE)
#' gr2 <- makeGRangesFromDataFrame(
#'     data.frame(chr = "chr1", start = 11:15, end = 11:15,
#'                strand = '*', mC = 1, uC = 1:5),
#'     keep.extra.columns = TRUE)
#'
#' answer <- poolFromGRlist(list(gr1, gr2), stat = 'sum', verbose = FALSE)
#'
#' @importFrom matrixStats rowMedians
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom S4Vectors mcols<-
#' @importFrom methods as
#'
#' @export
poolFromGRlist <- function(LR, stat="sum", num.cores=1, tasks=0L,
                           prob=FALSE, column=1L, verbose=TRUE, ...) {
   if (verbose)
       message("*** Building a unique GRanges object from the list...")
   if (inherits(LR, "list")) {
       LR <- try(as(LR, "GRangesList"))
   }
   if (inherits(LR, "GRangesList")) {
       if (prob) { ## Apply Fisher transformation
           message("* prob == TRUE. Applying Fisher transformation at column ",
                   column)
           LR = lapply(LR, function(GR) {
           mcols(GR) <- atanh(as.vector(mcols(GR)[, column]))
           return(GR)
           })
       }
       x0 <- uniqueGRanges(LR, num.cores=num.cores, tasks=tasks,
                       verbose=verbose, ...)
   } else {
       if (inherits(LR, "GRanges")) {
           x0 <- LR
       } else {
           stop(paste0("Object LR is neither a list of GRanges objects,",
                  " a 'GRangesList' object or a GRanges object"))
       }
   }

   if (verbose)  message("*** Building a virtual methylome...")
   x1 <- as.matrix(mcols(x0))
   cn <- colnames(mcols(x0))
   statist <- function(x, stat) {
       x <- switch(stat,
               sum=rowSums(x),
               mean=round(rowMeans(x)),
               median=round(rowMedians(x)))
   }

   if (prob) {
       ## Apply inverse of Fisher transformation
       prob <- statist(x1, stat="mean")
       mcols(x0) <- data.frame(prob=tanh(prob))
   } else {
       idx <- grep("mC", cn)
       mC <- statist(x1[, idx], stat=stat)
       idx <- grep("uC", cn)
       uC <- statist(x1[, idx], stat=stat)
       mcols(x0) <- data.frame(mC, uC)
   }
   return(x0)
}
