#' @rdname simulateCounts
#' @title Simulate read counts of methylated and unmethylated cytosines
#' @description Auxiliary function to simulate read counts of methylated and
#'     unmethylated cytosines
#'@details Methylation coverages (minimum 10) are generated from a Negative
#'     Binomial distribution with function \code{\link[MASS]{rnegbin}} from R
#'     package MASS. This function uses the representation of the Negative
#'     Binomial distribution as a continuous mixture of Poisson distributions
#'     with Gamma distributed means. Prior methylation levels are randomly
#'     generated with beta distribution using \code{\link[stats]{Beta}}
#'     function from R package “stats” and posterior methylation levels are
#'     generated according Bayes' theorem. The read of methylation counts are
#'     obtained as the product of coverage by the posterior methylation level.
#' @param num.samples Number of samples to generate.
#' @param sites Number of cytosine sites for each sample.
#' @param alpha Alpha parameter of beta distribution. Parameter shape1 from
#'     \code{\link[stats]{Beta}} function.
#' @param beta Beta parameter of beta distribution. Parameter shape2 from
#'     \code{\link[stats]{Beta}} function.
#' @param size number of trials (11 or more). Expected cytosine coverage.
#' @param theta Parameter theta from \code{\link[MASS]{rnegbin}}
#'     (overdispersion parameter).
#' @param sample.ids Names for the samples.
#' @importFrom stats rbeta rbinom
#' @importFrom MASS rnegbin
#' @return A list of GRanges objects with the methylated and unmethylated counts
#'     in its metacolumn.
#' @export
#' @author Robersy Sanchez
#' @examples
#' # *** Simulate samples with expected average of difference of methylation
#' # levels equal to 0.0427.
#' # === Expected mean of methylation levels ===
#' bmean <- function(alpha, beta) alpha/(alpha + beta)
#' bmean(0.03, 0.5) - bmean(0.007, 0.5) #' Expected difference = 0.04279707
#'
#' # === The number of cytosine sitesto generate ===
#' sites = 5000
#' # == Set a seed for pseudo-random number generation ===
#' set.seed(123)
#'
#' # === Simulate samples ===
#' ref = simulateCounts(num.samples = 1, sites = sites, alpha = 0.007,
#'                     beta = 0.5, size = 50, theta = 4.5, sample.ids = "C1")
#' treat = simulateCounts(num.samples = 2, sites = sites, alpha = 0.03,
#'                     beta = 0.5, size = 50, theta = 4.5,
#'                     sample.ids = c("T1", "T2"))
#'
#' #  === Estime Divergences ===
#' HD = estimateDivergence(ref = ref$C1, indiv =  treat, Bayesian = TRUE,
#'                         num.cores = 1L, percentile = 1)
#'
#' # === Difference of methylation levels of treatment simulated samples.
#' # Treatment versus reference
#' data.frame(mean.diff = c(mean(HD$T1$TV), mean(HD$T2$TV)),
#'            c("T1", "T2"), row.names = 2)

simulateCounts <- function(num.samples, sites, alpha, beta, size, theta,
                           sample.ids = NULL){
  coverage <- rnegbin(n = sites, mu = 100, theta = theta)
  # fix sites w/o reads
  coverage <- ifelse(coverage<10, 10, coverage)
  LR = list()
  for(k in 1:num.samples) {
    p = rbeta(n = sites, shape1 = alpha, shape2 = beta)
    shape1 <- theta * p
    shape2 <- theta * (1 - p)
    p = rbinom(n = sites, size=size,
               prob = rbeta(n=sites, shape1=shape1, shape2=shape2))/size
    mC = ceiling(coverage * p)
    uC = coverage - mC
    LR[[k]] = makeGRangesFromDataFrame(
      data.frame(chr = '1', start = 1:sites,
                 end = 1:sites, strand = '*',
                 mC = mC, uC = uC),
      keep.extra.columns = TRUE)
  }
  if (!is.null(sample.ids)) names(LR) <- sample.ids
  return(LR)
}
