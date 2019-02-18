#' @rdname getPotentialDIMP
#'
#' @title Potential methylation signal
#' @description This function perform a selection of the cytosine sites
#'     carrying the potential methylation signal. The potential signals from
#'     controls and treatments are used as prior classification in further step
#'     of signal detection.
#' @details The potential signals are cytosine sites k with information
#'     divergence (DIV_k) values greater than the DIV(alpha = 0.05). The value
#'     of alpha can be specified. For example, potential signals with
#'     DIV_k > DIV(alpha = 0.01) can be selected. For each sample, cytosine
#'     sites are selected based on the corresponding fitted Weilbull
#'     distribution model that has been supplied.
#'
#' @param LR A list of GRanges objects. Each GRanges object carrying information
#'     divergence values for each cytosine site in its meta-column.
#' @param nlms A list of distribution fitted models (output of
#'     'fitNonlinearWeibullDist' function) or NULL. If NULL, then empirical
#'     cumulative distribution function is used to get the potential DIMPs.
#' @param div.col Column number for divergence variable is located in the
#'     meta-column.
#' @param dist.name name of the distribution to fit: Weibull2P (default:
#'     "Weibull2P"), Weibull three-parameters (Weibull3P), gamma with
#'     three-parameter (Gamma3P), gamma with two-parameter (Gamma2P),
#'     generalized gamma with three-parameter ("GGamma3P") or four-parameter
#'     ("GGamma4P"), and the empirical cumulative distribution function (ECDF).
#' @param absolute Logic (default, FALSE). Total variation (TV, the difference
#'     of methylation levels) is normally an output in the downstream MethylIT
#'     analysis. If 'absolute = TRUE', then TV is transformed into |TV|, which
#'     is an information divergence that can be fitted to Weibull or to
#'     Generalized Gamma distribution. So, if the nonlinear fit was performed
#'     for |TV|, then absolute must be set to TRUE.
#' @param alpha A numerical value (usually alpha < 0.05) used to select
#'     cytosine sites k with information divergence (DIV_k) for which Weibull
#'     probability P[DIV_k > DIV(alpha)].
#' @param tv.col Column number for the total variation to be used for filtering
#'     cytosine positions (if provided).
#' @param tv.cut If tv.cut and tv.col are provided, then cytosine sites k with
#'     abs(TV_k) < tv.cut are removed before to perform the ROC analysis.
#' @param min.coverage Cytosine sites with coverage less than min.coverage are
#'     discarded. Default: 0
#'
#' @return A list of GRanges objects, each GRanges object carrying the selected
#'     cytosine sites and and the Weibull probability P[DIV_k > DIV(alpha)].
#'
#' @examples
#' num.points <- 1000
#' HD <- GRangesList( sample1 = makeGRangesFromDataFrame(
#'         data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
#'             strand = '*',
#'             hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
#'         keep.extra.columns = TRUE))
#' nlms <- nonlinearFitDist(HD, column = 1, verbose = FALSE)
#' getPotentialDIMP(LR = HD, nlms = nlms, div.col = 1, alpha = 0.05)
#'
#' @export
#'
getPotentialDIMP <- function(LR, nlms=NULL, div.col, dist.name = "Weibull2P",
                             absolute=FALSE, alpha=0.05, tv.col=NULL,
                             tv.cut=NULL, min.coverage=NULL) {

   P <- function(k) {
       d <- LR[[k]]

       if (!is.null(min.coverage)) {
         cov1 <- d$c1 + d$t1
         cov2 <- d$c2 + d$t2
         idx <- which((cov1 >= min.coverage) | (cov2 >= min.coverage))
         d <- d[idx]
       }
       q <- mcols(d[, div.col])[, 1]

       if (dist.name == "ECDF") ECDF <- ecdf(q)

       if (!is.null(tv.col) && !is.null(tv.cut)) {
           d <- d[ which( abs(mcols(d[, tv.col])[, 1]) > tv.cut)]
       }

       q <- mcols(d[, div.col])[, 1]
       if (absolute) q = abs(q)
       if (!is.null(nlms)) {
           m <- nlms[[k]]
           m <- m[, 1]
       } else {dist.name="ECDF"; ECDF <- ecdf(q)}

       if (dist.name != "ECDF") {
           p <- switch(dist.name,
                   Weibull2P=pweibull(q, shape=m[1], scale=m[2],
                                      lower.tail=FALSE),
                   Weibull3P=pweibull(q - m[3], shape=m[1], scale=m[2],
                                       lower.tail = FALSE),
                   Gamma2P=pgamma(q, shape=m[1], scale=m[2],
                                   lower.tail = FALSE),
                   Gamma3P=pgamma(q - m[3], shape=m[1], scale=m[2],
                                   lower.tail = FALSE),
                   GGamma3P=pggamma(q, alpha=m[1], scale=m[2], psi=m[3],
                                   lower.tail = FALSE),
                   GGamma4P=pggamma(q, alpha=m[1], scale=m[2], mu=m[3],
                                   psi=m[4], lower.tail = FALSE)
           )
       } else p <- (1 - ECDF(q))

       idx <- which(p < alpha)
       p <- p[idx]
       d <- d[idx]
       mcols(d) <- data.frame(mcols(d), wprob=p)
       return(d)
   }
   sn <- names(LR)
   LR <- lapply(1:length(LR), P)
   names(LR) <- sn
   return(LR)
}
