#' @rdname betaBinPosteriors
#'
#' @title Beta binomial posteriors
#' @description Here, beta binomial posteriors are estimated in the context of
#'     the methylation level estimation. The direct use of counts for
#'     methylation level estimation is not recommended, since there is external
#'     source of noise coming from the sample manipulation, the sequencing
#'     machine and sequence alignment that could alter the counts and the
#'     number of methylated cytosines with low or zero counts (mC = 0) is
#'     generally high. The low counts and noise imply that read counts of
#'     methylated (mC) and non-methylated (uC) cytosines must not be used
#'     directly in the estimation of methylation levels.
#'
#'     In a Bayesian framework, methylated read counts are modeled by a beta-
#'     binomial distribution, which accounts for both, the biological and
#'     sampling variations [1-3]. In our case we adopted the Bayesian approach
#'     suggested in reference [3](Chapter 3). Naive distribution q
#'     (methylation levels). In a Bayesian framework with uniform priors, the
#'     methylation level can be defined as:
#'     meth_level = ( mC + 1 )/( mC + uC + 2 ). However, the most
#'     natural statistical model for replicated BS-seq DNA methylation
#'     measurements is beta-binomial (the beta distribution is a prior
#'     conjugate of binomial distribution), we consider the p parameter
#'     (methylation level) in the binomial distribution as randomly drawn from
#'     a beta distribution. The hyper-parameters alpha ('a') and beta ('b')
#'     from the beta-binomial distribution are interpreted as pseudo-counts.
#' @details The posterior methylation levels are estimated as
#'     (a + success)/(a + b + trials), where 'a' and 'b' are the shape
#'     parameters of the beta distribution, alpha and beta parameters,
#'     respectively.
#'
#' @param success number of successful events
#' @param trials total number of events
#' @param a previous number of successful events
#' @param b previous number of unsuccessful events
#'
#' @return a probability
#'
#' @examples
#'     MethylIT:::.betaBinPosteriors(2, 8, 2, 8)
#'
#' @references 1. Hebestreit K, Dugas M, Klein H-U (2013) Detection of
#'     significantly differentially methylated regions in targeted bisulfite
#'     sequencing data. Bioinformatics 29: 1647-1653. Available:
#'     http://www.ncbi.nlm.nih.gov/pubmed/23658421. Accessed 4 February 2014.
#'     2. Robinson MD, Kahraman A, Law CW, Lindsay H, Nowicka M, et al. (2014)
#'     Statistical methods for detecting differentially methylated loci and
#'     regions. Front Genet 5: 324. Available:
#'     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4165320/ Accessed 26 August
#'     2015. 3. Baldi P, Brunak S (2001) Bioinformatics: the machine learning
#'     approach. Second. Cambridge: MIT Press. 452 p.
#' @keywords internal
.betaBinPosteriors <- function(success, trials, a, b) {
   (a + success)/(a + b + trials)
}
