#' @rdname estimateDivergence
#'
#' @title Information divergence estimator in respect to a reference sample
#' @description Wrapper of 'InfDiv' function to operate on list of GRanges
#' @details For the current version, the Information divergence of methylation
#'     levels is estimated based on Hellinger divergence (H). If read counts are
#'     provided, then Hellinger divergence is computed as given in the first
#'     formula from Theorem 1 from reference 1. In the present case:
#'
#'     \deqn{H = 2*(n[1] + 1)*(n[2] + 1)*((sqrt(p[1]) - sqrt(p[2]))^2 +
#'          (sqrt(1-p[1]) - sqrt(1-p[2]))^2)/(n[1] + n[2] + 2)}
#'
#'     where n[1] and n[2] are the coverage for the control and treatment,
#'     respectively. Notice that each row from the matrix of counts correspond
#'     to a single cytosine position and has four values corresponding to "mC1"
#'     and "uC1" (control), and mC2" and "uC2" for treatment.
#'
#'     According with the above equation, to estimate Hellinger divergence, not
#'     only the methylation levels are considered in the estimation of H,
#'     but also the control and treatment coverage at each given cytosine site.
#'     At this point, it is worthy to do mention that if the reference sample is
#'     derived with function \code{\link{poolFromGRlist}} using the 'sum' of
#'     read counts to conpute a methylation pool, then 'min.coverage' parameter
#'     value must be used to prevent an over estimation of the divergence for
#'     low coverage cytosines sites. For example, if a reference sample is
#'     derived as the methylation pool of read count sum from 3 individuals and
#'     we want to consider only methylation sites with minimum coverage of 4,
#'     then we can set min.coverage = c(12, 4), where the number 12 (3 x 4) is
#'     the minimum coverage requested for the each cytosine site in the
#'     reference sample.
#'
#'     If the methylation levels are provided in place of counts, then
#'     Hellinger divergence is computed as:
#'     \deqn{H = (sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1 - p[1]) -
#'           sqrt(1 - p[2]))^2}
#'
#'     This formula assumes that the probability vectors derived from the
#'     methylation levels (p_ij) p_j = c(p_ij, 1 - p_ij) (see function
#'     'estimateHellingerDiv') are an unbiased estimation of the expected one.
#'     The function applies a pairwise filtering after building a single GRanges
#'     from the two GRanges objects. Experimentally available cytosine sites are
#'     paired using the function 'uniqueGRanges'.
#'
#'     It is important to observe that several filtering conditions are provided
#'     to select biological meaningful cytosine positions, which prevent to
#'     carry experimental errors in the dowstream analyses. By filtering the
#'     read count we try to remove bad quality data, which would be in the edge
#'     of the experimental error originated by the BS-seq sequencing. It is
#'     responsability of the user to check whether cytosine positions used in
#'     the analysis are biological meaningful. For example, a cytosine position
#'     with counts mC1 = 10 and uC1 = 20 in the 'ref' sample and mC2 = 1 & uC2 =
#'     0 in an 'indv' sample will lead to methylation levels p1 = 0.333 and p2 =
#'     1, respectively, and TV = p2 - p1 = 0.667, which apparently indicates a
#'     hypermethylated site. However, there are not enough reads supporting p2 =
#'     1. A Bayesian estimation of TV will reveal that this site would be, in
#'     fact, hypomethylated. So, the best practice will be the removing of sites
#'     like that. This particular case is removed under the default settings:
#'     min.coverage = 4, min.meth = 4, and min.umeth = 0 (see example for
#'     function \code{\link{uniqueGRfilterByCov}}, called by
#'     estimateDivergence).
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
#' @param min.coverage An integer or an integer vector of length 2. Cytosine
#'     sites where the coverage in both samples, 'x' and 'y', are less than
#'     'min.coverage' are discarded. The cytosine site is preserved, however, if
#'     the coverage is greater than 'min.coverage'in at least one sample. If
#'     'min.coverage' is an integer vector, then the corresponding min coverage
#'     is applied to each sample.
#' @param min.meth An integer or an integer vector of length 2. Cytosine sites
#'     where the numbers of read counts of methylated cytosine in both samples,
#'     '1' and '2', are less than 'min.meth' are discarded. If 'min.meth' is an
#'     integer vector, then the corresponding min number of reads is applied to
#'     each sample. Default is min.meth = 4.
#' @param min.umeth An integer or an integer vector of length 2. Min number of
#'     reads to consider cytosine position. Specifically cytosine positions
#'     where (uC <= min.umeth) & (mC > 0) & (mC <= min.meth[1]) hold will be
#'     removed, where mC and uC stand for the numbers of methylated and
#'     unmethylated reads. Default is min.umeth = 0.
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
#' @return An object from "infDiv" class with the four columns of counts, the
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
estimateDivergence <- function(ref, indiv, Bayesian = FALSE, columns = NULL,
                           min.coverage = 4, min.meth = 4, min.umeth = 0,
                           high.coverage = NULL, percentile = 0.999,
                           num.cores = 1L, tasks = 0L,
                           meth.level = FALSE, verbose = TRUE, ...) {

   if (is.null(columns) && (!meth.level)) columns <- c(1,2)
   if (meth.level && (is.null(columns))) columns <- 1
   sn <- names(indiv)

   if (Sys.info()['sysname'] == "Linux") {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
   } else {
       bpparam <- SnowParam(workers = num.cores, type = "SOCK")
   }
   if (ncol(mcols(ref)) > 2) ref <- ref[ , columns]
   indiv <- lapply(indiv, function(x) x[, columns])

   if (meth.level) {
       x = bplapply(seq_len(length(indiv)), function(k, ref, indv, sn) {
           if (verbose) message("*** Processing sample #", k, " ", sn[k])
           x <- indv[[k]]
           x <- x[ ,columns]
           x <- uniqueGRanges(list(ref,x), num.cores=1L, tasks=tasks,
                               verbose=verbose, ...)
           x = estimateBayesianDivergence(x, Bayesian=Bayesian,
                               num.cores=1L, tasks=tasks,
                               meth.level=meth.level, verbose=verbose)
           return(x)
           }, BPPARAM=bpparam, ref=ref, indv=indiv, sn=sn)
   } else {
       x = bplapply(seq_len(length(indiv)), function(k, ref, indv, sn) {
           if (verbose) message("*** Processing sample #", k, " ", sn[ k ])
           x = uniqueGRfilterByCov(x=ref, y=indv[[k]],
                               min.coverage=min.coverage, min.meth = min.meth,
                               min.umeth = min.umeth, percentile=percentile,
                               num.cores=1L, tasks=tasks, verbose=verbose, ...)
           if (length(x) < 2)
               stop("*** At least two cytosine sites must pass the filtering ",
                   "conditions to estimate informations divergences. \n",
                   "The issue was found at sample number: ", k, ", id: ",
                   names(indv)[k])
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
