library(testthat)
library(GenomicRanges)
library(MethylIT)

context("nonlinearFitDist tests")

test_that("nonlinearFitDist with Weibull distribution values", {
  num.points <- 1000
  HD <- GRangesList(
    sample1 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
      keep.extra.columns = TRUE))

  nlms = nonlinearFitDist(HD, column = 1, verbose = FALSE)
  expect_true(as.numeric(as.character(nlms$sample1$R.Cross.val))[1] > 0.95)
})

test_that("nonlinearFitDist with non Weibull distribution values", {
  num.points <- 1000
  # Create data don't follow Hellinger divergence
  HD <- GRangesList(
    sample1 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = 0.75),
      keep.extra.columns = TRUE))

  #
  nlms = nonlinearFitDist(HD, column = 1, verbose = FALSE)
  # *** Trying nonlinear fit a Weibull 2P model...
  # *** Trying nonlinear fit a Weibull 3P model ...
  # Show Traceback
  #
  # Rerun with Debug
  # Error in `colnames<-`(`*tmp*`,
  #   value = c("Estimate", "Std. Error", "t value",  :
  #     attempt to set 'colnames' on an object with less than two dimensions
  #     In addition: Warning message:
  #       weibull3P(x, sample.size = sample.size,
  #                  npoints = npoints, npoints0 = npoints0,  :
  #       Data did not fit to the model.
  #       Returning empty coefficient table.
  TRUE
})

test_that("nonlinearFitDist with GGamma distribution values", {
  num.points <- 1000
  HD <- GRangesList(
    sample1 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rggamma(num.points, alpha = 0.41, psi = 1.6,
                                scale = 2.3)),
      keep.extra.columns = TRUE))
  nlms.w = nonlinearFitDist(HD, column = 1, verbose = FALSE)
  nlms.gg = nonlinearFitDist(HD, column = 1, dist.name = "GGamma3P",
                             verbose = FALSE)

  AIC.gg = as.numeric(as.character(nlms.gg$sample1$AIC)[1])
  AIC.w = as.numeric(as.character(nlms.w$sample1$AIC)[1])
  # To generate random values with GGamma disitribution.
  # Weilbull distribution is a particular case of GGamma. Therefore,
  # AIC.gg < AIC.w must hold.
  expect_true(AIC.gg < AIC.w)
})
