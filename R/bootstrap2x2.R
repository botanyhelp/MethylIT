#' @rdname bootstrap2x2
#'
#' @title bootstrap2x2
#' @description Parametric Bootstrap of 2x2 Contingence independence test. The
#'     goodness of fit statistic is the root-mean-square statistic (RMST) or
#'     Hellinger divergence, as proposed by Perkins et al. [1, 2]. Hellinger
#'     divergence (HD) is computed as proposed in [3].
#'
#' @param x A numerical matrix corresponding to cross tabulation (2x2) table
#'     (contingency table).
#' @param stat Statistic to be used in the testing: 'rmst','hdiv', or "all".
#' @param num.permut Number of permutations.
#'
#' @details For goodness-of-fit the following null hypothesis is tested
#'     \eqn{H_\theta: p = p(\theta)}
#'     To conduct a single simulation, we perform the following three-step
#'      procedure [1,2]:
#' \enumerate{
#'     \item To generate m i.i.d. draws according to the model distribution
#'           \eqn{p(\theta)}, where \eqn{\theta'} is the estimate calculated
#'           from the experimental data,
#'     \item To estimate the parameter \eqn{\theta} from the data generated in
#'           Step 1, obtaining a new estimate \eqn{\theta}est.
#'     \item To calculate the statistic under consideration (HD,
#'           RMST), using the data generated in Step 1 and taking the model
#'           distribution to be \eqn{\theta}est, where \eqn{\theta}est is the
#'           estimate calculated in Step 2 from the data generated in Step 1.
#' }
#'     After conducting many such simulations, the confidence level for
#'     rejecting the null hypothesis is the fraction of the statistics
#'     calculated in step 3 that are less than the statistic calculated from
#'     the empirical data. The significance level α is the same as a confidence
#'     level of \eqn{1-\alpha}.
#'
#' @return A p-value probability
#'
#' @examples
#'     set.seed(123)
#'     TeaTasting = matrix(c(8, 350, 2, 20), nrow = 2,
#'                         dimnames = list(Guess = c("Milk", "Tea"),
#'                         Truth = c("Milk", "Tea")))
#'     ## Small num.permut for test's speed sake
#'     bootstrap2x2( TeaTasting, stat = "all", num.permut = 100 )
#' @references
#' \enumerate{
#'     \item Perkins W, Tygert M, Ward R. Chi^2 and Classical Exact Tests
#'           Often Wildly Misreport Significance; the Remedy Lies in Computers
#'           [Internet]. Uploaded to ArXiv. 2011. Report No.:
#'           arXiv:1108.4126v2.
#'     \item Perkins, W., Tygert, M. & Ward, R. Computing the confidence
#'           levels or a root-mean square test of goodness-of-fit. 217,
#'           9072-9084 (2011).
#'     \item Basu, A., Mandal, A. & Pardo, L. Hypothesis testing for two
#'           discrete populations based on the Hellinger distance. Stat.
#'           Probab. Lett. 80, 206-214 (2010).
#' }
#' @export
bootstrap2x2 <- function(x, stat="rmst", num.permut=100) {
   HD <- function(x, y) {
       ## Function to compute the Hellinger divergence from
       ## an observed 2x2 contingence table
       v <- c(x,y)
       n1 <- sum(x)
       n2 <- sum(y)
       n <- cbind(n1, n2)
       ## pseudo-counts added if at least one of the cell counts is zero
       if (sum(v == 0) > 0) {
           p1 <- (x + 1) / (n1 + 2)
           p2 <- (y + 1) / (n2 + 2)
       } else {
           p1 <- x/n1
           p2 <- y/n2
       }
       p = cbind(p1, p2)
       hdiv <- (2 * (n[1] + 1) * (n[2] + 1) *
           sum((sqrt(p1) - sqrt(p2)) ^ 2 ) / (n[1] + n[2] + 2))
       return(hdiv)
   }

   m0 <- rowSums(x)
   n0 <- colSums(x)
   N0 <- sum(x)

   ## The expected number of counts
   freq0 <- as.vector(outer(m0, n0) / N0)
   prob <- freq0 / N0

   statis <- function(stat="rmst") {
       ## Function to randomly generate a 2x2 table based on the expected
       ## number of counts estimated from the observed 2x2 contingence table
       rtable <- rep(0, 4) ## all initial counts equal to zero
       ## Randomly select a cell in the table with probability equal to
       ## the expected counts divided by the total number of counts (N) in the
       ## table & increment the value in this cell by one. Repeat this
       ## procedure N times
       rcounts <- table(sample(x=1:4, size=N0, replace=TRUE,
                               prob=prob))
       ## updates initial counts
       rtable[ as.numeric(names(rcounts))] <- rcounts
       ## Compute the specified statistic for the randomly generated table
       y <- matrix(rtable, nrow=2)
       m <- rowSums(y)
       n <- colSums(y)
       N <- sum(y)
       freq <- as.vector(outer(m, n) / N)

       st <- switch(stat,
                   rmst=sum((rtable - freq) ^ 2) / 4,
                   hdiv=HD(rtable, freq),
                   all=list(rmst=sum((rtable - freq) ^ 2) / 4,
                               hdiv=HD(rtable, freq)))
       return(st)
   }

   x <- as.vector(x)
   if (stat == "all") {
       st0 <- c(rmst=sum((x - freq0) ^ 2) / 4, hdiv=HD(x, freq0))
       st <-  t(sapply(1:num.permut, function(i) statis(stat=stat)))
       res <- c(rmst.p.value=(sum(st[ ,1] > st0[1]) + 1) / num.permut,
               hdiv.p.value=(sum(st[ ,2] > st0[2]) + 1) / num.permut)
   } else {
       st0 <- switch(stat, rmst=sum((x - freq0) ^ 2) / 4, hdiv=HD(x, freq0))
       st <- sapply(1:num.permut, function(i) statis(stat=stat))
       res <- (sum(st > st0) + 1) / num.permut
   }
   return(res)
}
