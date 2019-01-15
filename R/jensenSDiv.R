#' @rdname jensenSDiv
#' @title Compute Jensen-Shannon Divergence
#' @description Compute Jensen-Shannon Divergence of probqbility vectors p
#'     and q.
#' @details The Jensen–Shannon divergence is a method of measuring the
#'     similarity between two probability distributions. Here, the
#'     generalization given in reference [1] is used. Jensen–Shannon divergence
#'     is expressed in terms of Shannon entroppy. 0 < jensenSDiv(p, q) < 1,
#'     provided that the base 2 logarithm is used in the estimation of the
#'     Shannon entropies involved.
#' @param p,q Probability vectors, sum(p_i) = 1 and sum(q_i) = 1.
#' @param Pi Weight of the probability distribution p. The weight for q is:
#'     1 - Pi. Default Pi = 0.5.
#' @param logbase A positive number: the base with respect to which logarithms
#      are computed (default: logbase = 2).
#' @examples
#' set.seed(123)
#' counts = sample.int(10)
#' prob.p = counts/sum(counts)
#' counts = sample.int(12,10)
#' prob.q = counts/sum(counts)
#' jensenSDiv(prob.p, prob.q)
#' @references
#' 1. J. Lin, “Divergence Measures Based on the Shannon Entropy,”
#' IEEE Trans. Inform. Theory, vol. 37, no. 1, pp. 145–151, 1991.
#' @export

jensenSDiv <- function(p, q, Pi=0.5, logbase=2) {
  m=Pi * p + (1 - Pi) * q
  Sm=shannonEntr(p=m, logbase=logbase)
  Sp=shannonEntr(p=p, logbase=logbase)
  Sq=shannonEntr(p=q, logbase=logbase)
  return(Sm - Pi * Sp - (1 - Pi) * Sq)
}
