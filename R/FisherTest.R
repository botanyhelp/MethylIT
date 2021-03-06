#' @rdname FisherTest
#' @title Fisher's exact test for read counts on GRanges objects
#' @description Given a GRanges object with the methylated and unmethylated
#'     read counts for control and treatment in its metacolumn, Fisher's exact
#'     test is performed for each cytosine site.
#' @details Samples from each group are pooled according to the statistic
#'     selected (see parameter pooling.stat) and a unique GRanges object is
#'     created with the methylated and unmethylated read counts for each group
#'     (control and treatment) in the metacolumn. So, a contingency table can be
#'     built for range from GRanges object.
#' @param LR A list of GRanges, a GRangesList, a CompressedGRangesList object,
#'     or an object from Methyl-IT downstream analyses: 'InfDiv' or "pDMP"
#'     object. Each GRanges object from the list must have two columns:
#'     methylated (mC) and unmethylated (uC) counts. The name of each element
#'     from the list must coincide with a control or a treatment name.
#' @param count.col 2d-vector of integers with the indexes of the read count
#'     columns. If not given, then it is assumed that the methylated and
#'     unmethylated read counts are located in columns 1 and 2 of each GRanges
#'     metacolumns. If object LR is the output of Methyl-IT function
#'     \code{\link[MethylIT]{estimateDivergence}}, then columns 1:4 are the read
#'     count columns: columns 1 and 2 are methylated and unmethylated read
#'     counts from the reference group, while columns 3 and 4 are methylated and
#'     unmethylated read counts from the treatment group, respectively. In this
#'     case, if the requested comparison is reference versus treatment, then no
#'     specification is needed for count.col. The comparison control versus
#'     treatment can be obtained by setting count.col = 3:4 and providing
#'     control.names and treatment.names.
#' @param control.names,treatment.names Names/IDs of control and treatment
#'     samples, which must be included in the variable GR at the metacolumn.
#'     Default NULL. If provided the Fisher's exact test control versus
#'     treatment is performed. Default is NULL. If NULL, then it is assumed that
#'     each GRanges object in LR has four columns of counts. The first two
#'     columns correspond to the methylated and unmethylated counts from
#'     control/reference and the other two columns are the methylated and
#'     unmethylated counts from treatment, respectively.
#' @param pooling.stat statistic used to estimate the methylation pool: row sum,
#'     row mean or row median of methylated and unmethylated read counts across
#'     individuals. If the number of control samples is greater than 2 and
#'     pooling.stat is not NULL, then they will pooled. The same for treatment.
#'     Otherwise, all the pairwise comparisons will be done.
#' @param tv.cut A cutoff for the total variation distance (TVD; absolute value
#'     of methylation levels differences) estimated at each site/range as the
#'     difference of the group means of methylation levels. If tv.cut is
#'     provided, then sites/ranges k with \eqn{|TV_k| < tv.cut} are removed
#'     before performing the regression analysis. Its value must be NULL or a
#'     number \eqn{0 < tv.cut < 1}.
#' @param hdiv.cut An optional cutoff for the Hellinger divergence (*hdiv*). If
#'     the LR object derives from the previous application of function
#'     \code{\link{estimateDivergence}}, then a column with the *hdiv* values is
#'     provided. If combined with tv.cut, this permits a more effective
#'     filtering of the signal from the noise. Default is NULL.
#' @param hdiv.col Optional. Columns where *hdiv* values are located in each
#'     GRanges object from LR. It must be provided if together with *hdiv.cut*.
#'     Default is NULL.
#' @param pAdjustMethod method used to adjust the results; default: BH
#' @param pvalCutOff cutoff used then a p-value adjustment is performed
#' @param saveAll if TRUE all the temporal results are returned
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see bpapply function from
#'     BiocParallel).
#' @param tasks integer(1). The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#' @param verbose if TRUE, prints the function log to stdout
#' @param progressbar logical(1). Enable progress bar
#' @param ... Additional parameters for function
#'     \code{\link[MethylIT]{uniqueGRanges}}.
#'
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom stats fisher.test
#' @return The input GRanges object with the columns of Fisher's exact test
#'     p-value, total variation (difference of methylation levels), and
#'     p-value adjustment.
#'
#' @examples
#' ## Get a dataset of Hellinger divergency of methylation levels
#' ## from the package
#' data(HD)
#' ### --- To get the read counts
#' hd <- lapply(HD, function(hd) {
#'    hd = hd[1:10,3:4]
#'    colnames(mcols(hd)) <- c("mC", "uC")
#'    return(hd)
#' })
#'
#' FisherTest(LR = hd,
#'        pooling.stat = NULL,
#'         control.names = "C1",
#'         treatment.names = "T1",
#'         pAdjustMethod="BH",
#'         pvalCutOff = 0.05,
#'         num.cores = 1L,
#'         verbose=FALSE)
#'
#' @seealso \code{\link[MethylIT.utils]{rmstGR}}
#' @export
FisherTest <- function(LR, count.col=c(1,2), control.names=NULL,
                       treatment.names=NULL, pooling.stat = "sum", tv.cut=NULL,
                       hdiv.cut=NULL, hdiv.col=NULL, pAdjustMethod="BH",
                       pvalCutOff=0.05, saveAll=FALSE, num.cores=1L, tasks=0L,
                       verbose=FALSE, progressbar = TRUE, ...) {

   if (any(!unlist(lapply(LR, function(GR) is(GR, "GRanges")))))
       stop("At least one element from 'LR' is not a 'GRanges' object")

   if (inherits(LR, c('InfDiv', "pDMP"))) {
       validateClass(LR)
   }

   if (!is.null(control.names)&&!is.null(treatment.names))
     LR = try(LR[c(control.names, treatment.names)], silent=TRUE)
   if (inherits(LR, "try-error"))
     stop("List's names does not match control & treatment names")

   # === Auxiliar function to perform FT ===
   ftest <- function(GR, num.cores = num.cores, tasks = tasks, ...) {
       count.matrix = as.matrix(mcols(GR))
       p1 <- count.matrix[, c(1,2)]
       p1 <- p1[, 1]/rowSums(p1)
       p2 <- count.matrix[, 3:4]
       p2 <- p2[, 1]/rowSums(p2)
       TV = p2 - p1; rm(p1, p2); gc()

       ind <- FALSE; idx <- c()
       if (!is.null(tv.cut)) idx.tv <- which(abs(TV) >= tv.cut)
       if (!is.null(hdiv.cut) && !is.null(hdiv.col))
           idx.hdiv <- which(mcols(GR[, hdiv.col])[,1] >= hdiv.cut)

       if (!is.null(tv.cut) && (!is.null(hdiv.cut) && !is.null(hdiv.col))) {
           idx <- unique(c(idx.tv, idx.hdiv))
           ind <- !is.na(idx)
       } else {
           if (!is.null(tv.cut)) {idx <- idx.tv; ind <- !is.na(idx)}
           if (!is.null(hdiv.cut) && !is.null(hdiv.col)) {
               idx <- idx.hdiv;
               ind <- !is.na(idx)
            }
       }

       if (sum(ind) > 0) {
           idx = idx[ind]
           count.matrix = as.matrix(mcols(GR[idx]))
       } else count.matrix = as.matrix(mcols(GR))

       if (nrow(count.matrix) > 1 ) {
           count.matrix <- count.matrix[, c(1,2,3,4)]
           sites <- nrow(count.matrix)
           GR$TV <- TV
           GR$pvalue <- rep(1, length(GR))
           GR$adj.pval <- rep(1, length(GR))
           count.matrix <- split(count.matrix, row(count.matrix))
           if (verbose)
             cat("*** Performing Fisher's exact test... \n
                 # of sites after filtering: ", sites, "\n")
           if (Sys.info()['sysname'] == "Linux") {
             bpparam <- MulticoreParam(workers=num.cores, tasks=tasks,
                                       progressbar = progressbar)
           } else {
             bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                   progressbar = progressbar)
           }
           pvals <- unname(unlist(bplapply(count.matrix, function(v) {
             fisher.test(matrix(as.integer(v), 2, byrow = TRUE))$p.value
           }, BPPARAM=bpparam)))
       } else {
           count.matrix = count.matrix[c(1,2,3,4)]
           pvals <- fisher.test(matrix(count.matrix, 2, byrow = TRUE))$p.value
       }

       if ((!is.null(tv.cut) | !is.null(hdiv.cut)) &&
           !saveAll && length(idx) > 0) {
           GR <- GR[idx]
           GR$pvalue <- pvals
           GR$adj.pval <- p.adjust(pvals, method=pAdjustMethod)
       } else {
           if ((!is.null(tv.cut) | !is.null(hdiv.cut)) &&
               saveAll && length(idx) > 0) {
               GR$pvalue[idx] <- pvals
               GR$adj.pval[idx] <- p.adjust(pvals, method=pAdjustMethod)
           }
           if (is.null(tv.cut) && is.null(hdiv.cut)) {
               GR$pvalue <- pvals
               GR$adj.pval <- p.adjust(pvals, method=pAdjustMethod)
           }
       }

       if (!is.null(pvalCutOff) && !saveAll) {
           GR <- GR[ GR$adj.pval < pvalCutOff ]
       }

       return(GR)
   }

   # This tests each sample against the reference
   if (is.null(control.names) || is.null(treatment.names)) {
       res <- lapply(LR, function(GR) ftest(GR, num.cores = num.cores,
                                            tasks = tasks, verbose = verbose) )
   }

   if (!is.null(control.names)&&!is.null(treatment.names)) {
       ctrl <- LR[control.names]

       treat <- LR[treatment.names]

       if (!is.null(pooling.stat)) {
           ctrl <- poolFromGRlist(ctrl, stat=pooling.stat,
                                   num.cores=num.cores, verbose=verbose)
           colnames(mcols(ctrl)) <- c("c1", "t1") # control counts

           treat <- poolFromGRlist(treat, stat=pooling.stat,
                                   num.cores=num.cores, verbose=verbose)
           colnames(mcols(treat)) <- c("c2", "t2") # treatment counts

           GR <- uniqueGRanges(list(ctrl, treat), verbose=verbose, ...)
           res <- ftest(GR, num.cores = num.cores, tasks = tasks,
                       verbose=verbose )
       } else {
           ctrl <- lapply(ctrl, function(GR) {
               GR <- GR[, count.col]
               colnames(mcols(GR)) <- c("c1", "t1") # Control counts
               return(GR)
           })

           treat <- lapply(treat, function(GR) {
               GR <- GR[, count.col]
               colnames(mcols(GR)) <- c("c2", "t2") # Treatment counts
               return(GR)
           })
           res = list()
           i = 1
           test.name = c()
           for(j in seq_len(length(ctrl))){
               for (k in seq_len(length(treat))) {
                   test.name = c(test.name, paste0(control.names[j], ".",
                                               treatment.names[k]))
                   GR <- uniqueGRanges(list(ctrl[[j]], treat[[k]]),
                                       verbose = verbose, ...)
                   if (verbose)
                       cat("*** Testing", paste0(control.names[k], " versus ",
                           treatment.names[j]), "\n")
                   res[[i]]=ftest(GR, num.cores = num.cores,
                                   tasks = tasks, verbose = verbose)
                   i = i + 1
               }
               names(res) <- test.name
           }
       }
   }
   if (!is.list(res)) res <- list(groupComparison = res)
   if (inherits(LR, "InfDiv") || inherits(LR, "pDMP")) cl <- class(LR)
   else cl <- class(res)
   res <- structure(res, class = c(cl, "testDMP"))
   return(res)
}
