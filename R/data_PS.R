#' Simulated dataset of potential DMPs used in examples
#'
#' Each individuals sample includes 10000 cytosine positions
#'
#' @format PS is an object from class "pDMP" carrying in the meta-columns the
#' following variables:
#' \describe{
#'     \item{p1}{methylation level from the reference sample}
#'     \item{p2}{methylation level from the treatment sample}
#'     \item{TV}{the total variation distance (difference of
#'               methylation levels)}
#'     \item{hdiv}{Hellinger divergence}
#'     \item{wprob}{the probabilities:
#'     \eqn{wprob = 1 - CDF probability}}
#' }
#'
#''PS' is an object from class "pDMP" carrying the same meta-columns as 'HD'
#'(dataset) plus the probabilities: \eqn{wprob = 1 - CDF probability}. 'PS'
#'was obtained with function \code{\link{getPotentialDIMP}}.
"PS"
