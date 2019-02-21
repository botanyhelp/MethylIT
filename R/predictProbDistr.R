#' @rdname predict.ProbDistr
#' @aliases predict.ProbDistrList
#' @title Predict function for probability distributions in Methyl-IT
#' @description This is an utility function to get predictions from the
#'     probability distributions models used in Methyl-IT: Weibull, Gamma, and
#'     generalized Gamma. Some times, after the nonlinear fit of any of the
#'     mentioned modelsm we would like to evaluate the model output.
#' @details Predictions are based on the best model fit returned by function
#'     \code{\link{nonlinearFitDist}}. The possible prediction are: *density*,
#'     *quantiles*, *random number* or *probabilities*.
#' @param nlm An object carrying the best nonlinear fit for a distribution model
#'     obtained with function \code{\link{nonlinearFitDist}}.
#' @param pred Type of prediction resquested: *density* ("dens"),*quantiles*
#'     ("quant"), *random number* ("rnum") or *probabilities* ("prob").
#' @param q numeric vector of quantiles, probabilities or an interger if
#'     pred = "rnum".
#' @param dist.name name of the distribution to fit: Weibull2P (default:
#'     "Weibull2P"), Weibull three-parameters (Weibull3P), gamma with
#'     three-parameter (Gamma3P), gamma with two-parameter (Gamma2P),
#'     generalized gamma with three-parameter ("GGamma3P") or four-parameter
#'     ("GGamma4P").
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @examples
#' set.seed(1)
#' num.points <- 1000
#' HD <- makeGRangesFromDataFrame(
#'   data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
#'             strand = '*',
#'             hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
#'   keep.extra.columns = TRUE)
#' nlms <- nonlinearFitDist(list(HD), column = 1, verbose = FALSE)
#'
#' x=seq(0.1, 10, 0.05)
#' y <- predict(nlms[[1]], pred="dens", q = x,
#'                 dist.name="Weibull2P")
#' y1 <- dweibull(x, shape = 0.75, scale = 1)
#' # The maximum difference between the "theoretical" and estimated densities
#' max(abs(round(y, 2) - round(y1, 2)))
#'
#' @importFrom stats dweibull pweibull qweibull rweibull dgamma pgamma qgamma
#'     rgamma
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @export
#'
predict.ProbDistr<- function(nlm, ...) UseMethod("predict", nlm)
predict.ProbDistr <- function(nlm, pred="quant", q=0.95, dist.name) {
   if (missing(nlm)) stop("A probability distribution model must be provided")
   m <- nlm[, 1]

   arg.num <- switch(dist.name,
                       Weibull2P = length(m) == 2,
                       Weibull3P = length(m) == 3,
                       Gamma2P = length(m) == 2,
                       Gamma3P = length(m) == 3,
                       GGamma3P = length(m) == 3,
                       GGamma4P = length(m) == 4
                     )

   if (!arg.num) stop("The number of model parameters does not match the model")

   if (pred == "dens") {
       res <- try(switch(dist.name,
                       Weibull2P=dweibull(q, shape=m[1], scale=m[2]),
                       Weibull3P=dweibull(q, shape=m[1], scale=m[2]) + m[3],
                       Gamma2P=dgamma(q, shape=m[1], scale=m[2]),
                       Gamma3P=dgamma(q, shape=m[1], scale=m[2]) + m[3],
                       GGamma3P=dggamma(q, alpha=m[1], scale=m[2], psi=m[3]),
                       GGamma4P=dggamma(q, alpha=m[1], scale=m[2], mu=m[3],
                                       psi=m[3])
       ), silent = TRUE)
   }

   if (pred == "quant") {
       res <- try(switch(dist.name,
                       Weibull2P=qweibull(q, shape=m[1], scale=m[2]),
                       Weibull3P=qweibull(q, shape=m[1], scale=m[2]) + m[3],
                       Gamma2P=qgamma(q, shape=m[1], scale=m[2]),
                       Gamma3P=qgamma(q, shape=m[1], scale=m[2]) + m[3],
                       GGamma3P=qggamma(q, alpha=m[1], scale=m[2], psi=m[3]),
                      GGamma4P=qggamma(q, alpha=m[1], scale=m[2], mu=m[3],
                                       psi=m[3])
       ), silent = TRUE)
   }

   if (pred == "prob") {
       res <- try(switch(dist.name,
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
       ), silent = TRUE)
   }

   if (pred == "rnum") {
       res <- try(switch(dist.name,
                       Weibull2P=rweibull(q, shape=m[1], scale=m[2]),
                       Weibull3P=rweibull(q, shape=m[1], scale=m[2]) + m[3],
                       Gamma2P=rgamma(q, shape=m[1], scale=m[2]),
                       Gamma3P=rgamma(q, shape=m[1], scale=m[2]) + m[3],
                       GGamma3P=rggamma(q, alpha=m[1], scale=m[2], psi=m[3]),
                       GGamma4P=rggamma(q, alpha=m[1], scale=m[2], mu=m[3],
                                       psi=m[4])
     ), silent = TRUE)
   }

   if (!inherits(res, "try-error")) return(res) else {
       warning("Your model paramters return arror")
       return(res <- NA)
   }
}

#' @name predict.ProbDistrList
#' @rdname predict.ProbDistr
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @export

predict.ProbDistrList<- function(nlm, ...) UseMethod("predict", nlm)
predict.ProbDistrList <- function(nlm, pred="quant", q=0.95, dist.name,
                                  num.cores= 1L, tasks=0L) {
   if (Sys.info()['sysname'] == "Linux") {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
   } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")

   res <- bplapply(nlm, function(model) {
               predict(model, pred=pred, q=q, dist.name=dist.name)
           }, BPPARAM=bpparam
   )
   names(res) <- names(nlm)
   return(res)
}







