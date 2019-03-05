#' @rdname estimateDivergence
#'
#' @title Information divergence estimator in respect to a reference sample
#' @description Wrapper of 'InfDiv' function to operate on list of GRanges
#' @details For the current version, the Information divergence of methylation
#'     levels is estimated based on Hellinger divergence. If read counts are
#'     provided, then Hellinger divergence is computed as given in the first
#'     formula from Theorem 1 from reference 1. In the present case: hdiv =
#'     2*(n[1] + 1)*(n[2] + 1)*((sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1-p[1]) -
#'     sqrt(1-p[2]))^2)/(n[1] + n[2] + 2)
#'
#'     where n[1] and n[2] are the coverage for the control and treatment,
#'     respectively. Notice that each row from the matrix of counts correspond
#'     to a single cytosine position and has four values corresponding to "mC1"
#'     and "uC1" (control), and mC2" and "uC2" for treatment.
#'
#'     If the methylation levels are provided in place of counts, then
#'     Hellinger divergence is computed as:
#'     hdiv = (sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1 - p[1]) - sqrt(1 - p[2]))^2
#'
#'     This formula assumes that the probability vectors derived from the
#'     methylation levels (p_ij) p_j = c(p_ij, 1 - p_ij) (see function
#'     'estimateHellingerDiv') are an unbiased estimation of the expected one.
#'     The function applies a pairwise filtering after building a single GRanges
#'     from the two GRanges objects. Experimentally available cytosine sites are
#'     paired using the function 'uniqueGRanges'.
#'
#' @param ref The GRanges object of the reference individual that will be used
#'     in the estimation of the information divergence.
#' @param indiv A list of GRanges objects from the individuals that will be
#'     used in the estimation of the information divergence.
#'@param Bayesian Logical. Whether to perform the estimations based on
#'     posterior estimations of methylation levels.
#' @param columns Vector of one or two integer numbers denoting the indexes of
#'     the columns where the methylated and unmethylated read counts are found
#'     or, if meth.level = TRUE, the columns corresponding to the methylation
#'     levels. If columns = NULL and meth.level = FALSE, then columns = c(1,2)
#'     is assumed. If columns = NULL and meth.level = TRUE, then columns = 1 is
#'     assumed.
#' @param min.coverage Cytosine sites with coverage less than min.coverage are
#'     discarded.
#' @param high.coverage An integer for read counts. Cytosine sites having
#'     higher coverage than this are discarded.
#'@param percentile Threshold to remove the outliers from each file and all
#'     files stacked.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see 'bplapply' function from
#'     BiocParallel package).
#' @param tasks integer(1). The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#'@param meth.level Logic. Whether methylation levels are given in place of
#'     counts.
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... Additional parameters for 'uniqueGRanges' function.
#'
#' @return A list of GRanges objects with the four columns of counts, the
#'     information divergence, and additional columns: 1) The original matrix
#'     of methylated (c_i) and unmathylated (t_i) read counts from control
#'     (i=1) and treatment (i=2) samples. 2) p1" and "p2": methylation levels
#'     for control and treatment, respectively. 3) "bay.TV": total variation
#'     TV = p2 - p1. 4) "TV": total variation based on simple counts:
#'     TV=c1/(c1+t1)-c2/(c2+t2). 5) "hdiv": Hellinger divergence. If
#'     Bayesian = TRUE, the results are based on the posterior estimations of
#'     methylation levels.
#' @author Robersy Sanchez
#' @examples
#'     num.samples <- 250
#'     x <- data.frame(chr = "chr1", start = 1:num.samples,
#'                     end = 1:num.samples,strand = '*',
#'                     mC = rnbinom(size = num.samples, mu = 4, n = 500),
#'                     uC = rnbinom(size = num.samples, mu = 4, n = 500))
#'     y <- data.frame(chr = "chr1", start = 1:num.samples,
#'                     end = 1:num.samples, strand = '*',
#'                     mC = rnbinom(size = num.samples, mu = 4, n = 500),
#'                     uC = rnbinom(size = num.samples, mu = 4, n = 500))
#'     x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
#'     y <- makeGRangesFromDataFrame(y, keep.extra.columns = TRUE)
#'     HD <- estimateDivergence(ref = x, indiv = list(y))
#'
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom GenomicRanges GRanges GRangesList
#'
#' @export
estimateDivergence <- function(ref, indiv, Bayesian=FALSE, columns=NULL,
                          min.coverage=4, high.coverage=NULL, percentile=0.999,
                          num.cores=1L, tasks=0L, meth.level=FALSE,
                          verbose=TRUE, ...) {

   if (is.null(columns) && (!meth.level)) columns <- 1:2
   if (meth.level && (is.null(columns))) columns <- 1
   sn <- names(indiv)

   if (Sys.info()['sysname'] == "Linux") {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
   } else {
       bpparam <- SnowParam(workers = num.cores, type = "SOCK")
   }
   if (meth.level) {
       x = bplapply(1:length(indiv), function(k, ref, indv, sn) {
           if (verbose) message("*** Processing sample #", k, " ", sn[k])
           x <- indv[[k]]
           x <- x[ ,columns]
           ref <- ref[ ,columns]
           x <- uniqueGRanges(list(ref,x), num.cores=1L, tasks=tasks,
                               verbose=verbose, ...)
           x = estimateBayesianDivergence(x, Bayesian=Bayesian,
                               num.cores=1L, tasks=tasks,
                               meth.level=meth.level, verbose=verbose)
           return(x)
           }, BPPARAM=bpparam, ref=ref, indv=indiv, sn=sn)
   } else {
       x = bplapply(1:length(indiv), function(k, ref, indv, sn) {
           if (verbose) message("*** Processing sample #", k, " ", sn[ k ])
           x = uniqueGRfilterByCov(x=ref, y=indv[[k]],
                               min.coverage=min.coverage,
                               percentile=percentile, columns=columns,
                               num.cores=1L, tasks=tasks,
                               verbose=verbose, ...)
           x = estimateBayesianDivergence(x, Bayesian=Bayesian, num.cores=1L,
                               tasks=tasks, meth.level=meth.level,
                               verbose=verbose)
           return(x)
           }, BPPARAM=bpparam, ref=ref, indv=indiv, sn=sn)
   }
   names(x) <- sn
   x <- structure(x, class=c("InfDiv", "list"))
   return(x)
}
