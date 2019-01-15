#' @rdname ksTest
#'
#' @title Kolmogorov-Smirnov statistics
#' @description Permutation test for Kolmogorov-Smirnov statistics
#'
#' @param x numerical vector to perform the goodness of fit
#' @param CDF the name of the cummulative distribution function (CDF)
#' @param pars vector of parameters to evaluate the CDF:
#'     4P GG distribution: c(shape=value, scale=value, mu=value, psi=value)
#'     3P GG distribution: c(shape=value, scale=value, psi=value)
#'     3P Weibull distribution: c(shape=value, scale=value, mu=value)
#'     2P Weibull distribution: c(shape=value, scale=value)
#' @param num.sampl number of elements to be sampled
#' @param sample.size number of permutations. If sample.size < length(x), then
#'     the test becomes a Monte Carlo test
#' @param numcores number of cores
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... other parameters
#'
#' @return gamma distribution CDF
#'
#' @author Robersy Sanchez - 02/29/2016
#' @references Alastair Sanderson. Using R to analyse data statistical and
#'     numerical data analysis with R
#'     http://www.sr.bham.ac.uk/~ajrs/R/r-analyse_data.html
#'
#' @examples
#' num.samples <- 1000
#' x <- rweibull(num.samples, shape = 1.01, scale = 1.01)
#' ksTest(x, pars = c(shape = 1, scale = 1))
#'
#' @importFrom stats ks.test na.omit
#' @importFrom parallel mclapply
#' @export
ksTest <- function(x, CDF="Weibull", pars, num.sampl=999, sample.size,
                   numcores=1, verbose=TRUE, ... ) {

   funLIST <- c("pweibull", "pweibull", "pggamma3P", "pggamma3P")
   distNAMES <- c("Weibull", "Weibull 3P", "GGamma 3P", "GGamma 4P")
   ind <- as.numeric(na.omit( match(CDF, distNAMES)))
   distNAME <- funLIST[ind]

   ## All the estimations are based on gamma and Weibull distritions.
   ## So, a tranformation is applied to variable x.
   if (CDF == "Weibull 3P") {
       x <- x - pars[3]
       pars <- pars[1:2]
   }
   if (CDF == "GGamma 3P") {
       x <- (x/pars[2])^pars[1]
       pars <- pars[3]
   }
   if (CDF == "GGamma 4P") {
       x <- ((x - pars[3]) / pars[2]) ^ pars[1]
       pars <- pars[4]
   }

   if (missing(sample.size) || sample.size >= round(length(x)))
       sample.size <- round(length(x) / 3)

   kstest <- function(x, R, szise, mc.cores, distNAME) {
       myfun <- function(a) {
           suppressWarnings(unname(ks.test(a, distNAME, pars)$statistic))
       }

       DoIt <- function(r) {
           i <- sample(length(x), szise)
           myfun(x[i])  ## to test empirical versus theoretical values
       }

       pstats <- mclapply(1:R, DoIt, mc.cores=mc.cores)
       pstats <- unlist(pstats)
       stat <- myfun(x)
       list(p.value=mean(c(stat, pstats) >= stat, na.rm=TRUE),
           KS.stat=stat, boot.ks=pstats)
   }

   if (verbose) message( "*** Monte Carlo KS test,..\n" )
   kstest(x, R=num.sampl, szise=sample.size,
         mc.cores=numcores, distNAME=distNAME)
}
