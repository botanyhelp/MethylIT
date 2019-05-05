#' Simulated dataset used in the examples
#'
#' Each individuals sample includes 5000 cytosine position
#'
#' @format HD is an object from class "InfDiv" carrying:
#'     \describe{
#'         \item{p1}{methylation level from the reference sample}
#'         \item{p2}{methylation level from the treatment sample}
#'         \item{TV}{the total variation distance (difference of
#'                 methylation levels)}
#'         \item{hdiv}{Hellinger divergence}
#'     }
#'
#'     HD was obtained with function \code{\link{estimateDivergence}}
#' @format  nlms carries the information on the best fitted probability
#'     distribution model for each individual sample. 'nlms' was obtianed with
#'     function \code{\link{nonlinearFitDist}}.
#' @format PS is an object from class "pDMP" carrying the same meta-columns as
#'     HD plus the probilities: \eqn{wprob = 1 - Weibull probability}. 'PS' was
#'     obtained with function \code{\link{estimateCutPoint}}.
#'
#'
