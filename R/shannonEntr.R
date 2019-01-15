#' @rdname shannonEntr
#' @title Compute Shannon Entropy
#' @description Compute Shannon Entropy of probability vector p.
#' @details By definition, if p_i = 0 for some i, the value of the corresponding
#'     summ and 0*log(0) is taken to be 0.
#' @param p A probability vector, sum(p) = 1.
#' @param logbase A positive number: the base with respect to which logarithms
#      are computed (default: logbase = 2).
#' @examples
#' counts = sample.int(10)
#' prob = counts/sum(counts)
#' shannonEntr(prob)
#' @export
shannonEntr <- function(p, logbase = 2) {
   logb <- function(p) {
       n <- length(p)
       if (n > 1) {
           logP <- integer(n)
           idx <- p > 0
           logP[idx] <- -log(p[idx], base=logbase)
       } else logP <- -log(p, base=logbase)
       return(logP)
   }
   return(p * logb(p) + (1 - p) * logb(1 - p))
}


