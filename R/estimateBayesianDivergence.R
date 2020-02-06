#' @rdname estimateBayesianDivergence
#'
#' @title Information divergence estimator
#' @description The Information divergence of methylation levels is estimated
#'     using the direct estimation or a Bayesian approach of the methylation
#'     levels. Hellinger divergence is computed as given in reference 1.
#' @details For the current version, the Information divergence of methylation
#'     levels is estimated based on Hellinger divergence. If read counts are
#'     provided, then Hellinger divergence is computed as given in the first
#'     formula from Theorem 1 from reference 1. In the present case,
#'
#'     hdiv = 2*(n[1] + 1)*(n[2] + 1)*((sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1 -
#'     p[1]) - sqrt(1 - p[2]))^2)/(n[1] + n[2] + 2)
#'
#'     where n[1] and n[2] are the coverage for the control and treatment,
#'     respectively. Notice that each row from the matrix of counts correspond
#'     to a single cytosine position and has four values corresponding to "mC1"
#'     and "uC1" (control), and mC2" and "uC2" for treatment.
#'
#'     If the methylation levels are provided in place of counts, then
#'     Hellinger divergence is computed as:
#'     hdiv = (sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1 -p[1]) - sqrt(1 - p[2]))^2
#'
#'     This formula assumes that the probability vectors derived from the
#'     methylation levels (p_ij) p_j = c(p_ij, 1 - p_ij) (see function
#'     'estimateHellingerDiv') are an unbiased estimation of the expected one.
#'
#' @param x A matrix of counts or GRanges object with the table of counts in
#'     the meta-columns (methylated mC and unmethylated uC cytosines). Unless
#'     specified in the parameter 'columns', the methylation counts must be
#'     given in the first four columns: "mC1" and "uC1" methylated and
#'     unmethylated counts for control sample, and "mC2" and "uC2" methylated
#'     and unmethylated counts for treatment sample, respectively.
#' @param Bayesian logical. Whether to perform the estimations based on
#'     posterior estimations of methylation levels.
#' @param JD Logic (Default:FALSE). Option on whether to add a column with
#'     values of J-information divergence (see \code{\link{estimateJDiv}}).
#'     It is only compute if JD = TRUE and meth.level = FALSE.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see 'bplapply' function from
#'     BiocParallel package).
#' @param tasks integer(1). The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#' @param columns Vector of integer numbers of the columns where the counts
#'     "mC1", "uC1", "mC2", and "uC2" are in the matrix (default 1 to 4). That
#'     is, the input could have more than 4 columns, but only 4 columns with
#'     the counts are used.
#' @param meth.level methylation levels can be provided in place of counts.
#' @param preserve.gr Logic (Default:FALSE). Option of whether to preserve all
#'     the metadata from the original GRanges object.
#' @param logbase Logarithm base used to compute the JD (if JD = TRUE).
#'     Logarithm base 2 is used as default (bit unit). Use logbase = exp(1) for
#'     natural logarithm.
#' @param verbose if TRUE, prints the function log to stdout
#'
#' @return The input matrix or GRanges object with the four columns of counts
#'     and additional columns. If Bayessian = TRUE, the results are based on
#'     the posterior estimations of methylation levels. 1) p1" and "p2":
#'     methylation levels for control and treatment; 2) "hdiv": Hellinger
#'     divergence; 3) "bay.TV". "TV": total variation TV = p2 - p1 is based
#'     on simple counts
#'
#' @author Robersy Sanchez
#' @examples
#' ## The read count data are created
#'     x <- data.frame(chr = "chr1", start = 1:10,
#'                     end = 1:10,strand = '*',
#'                     mC1 = rnbinom(size = 10, mu = 4, n = 500),
#'                     uC1 = rnbinom(size = 10, mu = 4, n = 500),
#'                     mC2 = rnbinom(size = 10, mu = 4, n = 500),
#'                     uC2 = rnbinom(size = 10, mu = 4, n = 500))
#'     x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
#' ## Estimation of the information divergences
#'     hd <- estimateBayesianDivergence(x, JD = TRUE)
#'
#' ## Keep in mind that Hellinger and J divergences are, in general, correlated!
#'     cor.test(x = as.numeric(hd$hdiv), y = as.numeric(hd$jdiv),
#'             method = "kendall")
#'
#' @references 1. Basu  A., Mandal  A., Pardo L (2010) Hypothesis testing for
#'     two discrete populations based on the Hellinger distance. Stat Probab
#'     Lett 80: 206-214.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom S4Vectors mcols<-
#' @export
#' @keywords internal
estimateBayesianDivergence <- function(x, Bayesian=FALSE, JD = FALSE,
                                       num.cores=1, tasks=0L,
                                       columns=c(mC1=1, uC1=2, mC2=3, uC2=4),
                                       meth.level=FALSE,
                                       preserve.gr = FALSE,
                                       logbase = 2,
                                       verbose=TRUE) {

   if (Sys.info()['sysname'] == "Linux") {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
   } else {
       bpparam <- SnowParam(workers = num.cores, type = "SOCK")
   }
   ismatrix <- TRUE
   if (class(x) != "matrix") {
       HDiv <- x
       x <- as.matrix(mcols(x))
       ismatrix <- FALSE
   }
   ## If x carriers a read counts, then:
   if (!meth.level) {
       x <- x[ ,columns]
       r0 <- rowSums(x)
       ind <- which(r0 > 4); rm(r0)
       x <- x[ind, ]
       if (!ismatrix) HDiv <- HDiv[ind, ]
       n1 <- x[ ,1] + x[ ,2]
       n2 <- x[ ,3] + x[ ,4]
       n <- cbind(n1, n2)
       x0 <-  x
       p1 <- x[ ,1] / n1
       p2 <- x[ ,3] / n2

       ## if the coverage is zero in control or the reference
       ## individual, then p is NaN (NA). By definition these
       ## cytosine sites has methylation levels p = 0.
       p1[is.na(p1)] <- 0
       p2[is.na(p2)] <- 0

       ## Ordinary TV
       TV <- p2 - p1

       if (Bayesian) {
           if (nrow(x) < 10)
               stop(paste("*** You must provide at least 10 cytosine sites ",
                           "to apply a Bayessian approach \n",
                           "using beta distributed priors"))
           if (verbose) cat( "*** Estimating betaBinomial-posteriors... \n" )
           ## Naive distribution q (methylation levels).
           ## In a Bayesian framework with uniform priors,
           ## the methylation level can be defined as:
           ## meth_level = ( mC + 1 )/( mC + uC + 2 ).
           q1 <- (x[ ,1] + 1) / (n1 + 2)
           q2 <- (x[ ,3] + 1) / (n2 + 2)

           ## The shape parameters estimated with 'nlm'
           beta1 <- .estimateBetaDist(q1)
           beta2 <- .estimateBetaDist(q2)
           ## Assuming beta priors
           n1[n1 == 0] <- 2
           n2[n2 == 0] <- 2

           p1 <- .betaBinPosteriors(x[ ,1], n1, a=beta1[1], b=beta1[2])
           p2 <- .betaBinPosteriors(x[ ,3], n2, a=beta2[1], b=beta2[2])
       }

       x <- cbind(p1, p2)

       if (verbose) cat("*** Estimating Hellinger divergence... \n")
       hdiv <- bplapply(seq_len(nrow(x)), function(i) {
                                   estimateHellingerDiv(p=as.numeric(x[i, ]),
                                                        n=as.numeric(n[i, ]))},
                        BPPARAM=bpparam)
       if (verbose) cat( "* Coercing from list to vector...\n" )
       hdiv <- unlist(hdiv)
       if (Bayesian) {
           x <- data.frame(x0, p1, p2, TV, p2 - p1, hdiv)
           colnames(x) <- c("c1", "t1", "c2", "t2", "p1", "p2", "TV", "bay.TV",
                               "hdiv")
       } else {
           x <- data.frame(x0, p1, p2, TV, hdiv)
           colnames(x) <- c("c1", "t1", "c2", "t2", "p1", "p2", "TV", "hdiv")
       }
       if (JD) {
           jdiv <- bplapply(seq_len(nrow(x)), function(i) {
                           estimateJDiv(p = as.numeric(x[i, c("p1", "p2")]),
                                       logbase = logbase)},
                           BPPARAM=bpparam)
           x$jdiv <- unlist(jdiv)
       }
   } else {
       if (verbose) cat("*** Estimating Hellinger divergence... \n")
       hdiv <- bplapply(seq_len(nrow(x)), function(i) {
           estimateHellingerDiv(p=as.numeric(x[i, ]))}, BPPARAM=bpparam)
       if (verbose) cat( "* Coercing from list to vector...\n" )
       hdiv <- unlist(hdiv)
       x <- data.frame(x, TV=x[ ,2] - x[ ,1], hdiv)
       colnames(x) <- c("p1", "p2", "TV", "hdiv")
   }
   if (!ismatrix) {
       if (preserve.gr) {
           if (JD && !meth.level)
               mcols(HDiv) <- data.frame(mcols(HDiv),
                                       x[, c("p1", "p2", "TV", "hdiv", "jdiv")])
           else mcols(HDiv) <- data.frame(mcols(HDiv),
                                           x[, c("p1", "p2", "TV", "hdiv")])
       } else mcols(HDiv) <- x
       return(HDiv)
   } else return(x)
}
