#' @rdname estimateJDiv
#'
#' @title J Information Divergence of Methylation Levels
#' @description Given a the methylation levels from two individuals at a given
#'     cytosine site, this function computes the J information divergence
#'     (\emph{JD}) between methylation levels. The motivation to introduce
#'     \emph{JD} in Methyl-IT is founded on:
#'     \enumerate{
#'         \item It is a symmetrised form of Kullback–Leibler divergence
#'             (\eqn{D_KL}). Kullback and Leibler themselves actually defined
#'             the divergence as: \eqn{D_KL(P || Q) + D_KL(Q || P)}, which is
#'             symmetric and nonnegative, where the probability distributions
#'             \emph{P} and \emph{Q} are defined on the same probability space
#'             (see
#'             \href{https://is.gd/oS8Bwv}{Wikipedia}).
#'         \item In general, \emph{JD} is highly correlated with Hellinger
#'             divergence, which is the main divergence currently used in
#'             Methyl-IT (see examples for function
#'             \code{\link{estimateDivergence}}.
#'         \item By construction, the unit of measurement of \emph{JD} is:
#'             bit of information, which set the basis for further
#'             information-thermodynamics analyses.
#'     }
#'
#' @details The methylation level \eqn{p_ij} at a given cytosine site \emph{i}
#'     from an individual \emph{j} corresponds to a probability vector \eqn{q_ij
#'     = c(p_ij, 1 - p_ij)}. Then, the J-information divergence between the
#'     methylation levels \emph{p_i1} and \emph{p_i2} is the divergence between
#'     the vectors \eqn{q_i1 = (p_i1, 1 - p_i1)} and \eqn{q_i2 = c(p_i2, 1 -
#'     p_i2)}. \emph{JD} is computed as in reference (1):
#'
#' \deqn{JD(q_i1, q_i2) = (q_i1[1] * log(q_i1[1]/q_i2[1]) + q_i1[2] *
#'                         log(q_i1[2]/q_i2[2]) +
#'                 q_i2[1] * log(q_i2[1]/q_i1[1]) + q_i2[2] *
#'                             log(q_i2[2]/q_i1[2]))/2}.
#'
#' @param p A numerical vector of the methylation levels p = c(p1, p2) from
#'     individuals 1 and 2.
#' @param logbase Logarithm base used to compute the JD. Logarithm base 2 is
#'     used as default. Use logbase = exp(1) for natural logarithm.
#' @return The J divergence value for the given methylation levels is
#'     returned
#' @export
#' @examples
#' p <- c(0.5, 0.5)
#' estimateJDiv(p)
#'
#' ## A numerical trick is implemented. The J-divergence values are the
#' ## same for the following vectors:
#' p <- c(0.9999999999999999, 0)
#' q <- c(1, 0)
#'
#' estimateJDiv(p) == estimateJDiv(q)
#'
#' @author Robersy Sanchez 11/27/2019 <https://github.com/genomaths>
#' @seealso \url{https://en.wikipedia.org/wiki/Kullback-Leibler_divergence}
#'     for more details, and \code{\link{estimateDivergence}} for an example of
#'     using it.
#'
#' @references
#' \enumerate{
#'     \item Lin J. Divergence Measures Based on the Shannon Entropy. IEEE
#'           Trans Inform Theory, 1991, 37:145–51.
#'     \item Sanchez R, Mackenzie SA. Information thermodynamics of cytosine
#'           DNA methylation. PLoS One, 2016, 11:e0150427.
#' }
#'
#' @export
estimateJDiv <- function(p, logbase = 2) {
   if (any(p > 1) | any(p < 0))
       stop("*** Vector p has values out of the range [0, 1]")
   jdiv <- 0
   if (!is.na(sum(p)) && (sum(p) > 0)) {
      p[ p == 1] <- 0.9999999999999999
      p1 <- c(p[1], 1 - p[1])
      p2 <- c(p[2], 1 - p[2])

      jdiv <- (p1[1] * .log(p1[1]/p2[1], logbase = logbase) +
                   p1[2] * .log(p1[2]/p2[2], logbase = logbase) +
                   p2[1] * .log(p2[1]/p1[1], logbase = logbase) +
                   p2[2] * .log(p2[2]/p1[2], logbase = logbase))/2
   }
   return(jdiv)
}


## --- Auxiliary function  ------
.log <- function(p, logbase = 2) {
   logb <- function(p) {
      n <- length(p)
      logP <- integer(n)
      idx <- (p > 0 & p != Inf)
      logP[idx] <- log(p[idx], base = logbase)
      return(logP)
   }
   return(logb(p))
}

