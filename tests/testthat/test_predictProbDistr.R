library(testthat)
library(GenomicRanges)
library(MethylIT)

context("predictProbDistr tests")

test_that("predictProbDistr test", {
  set.seed(1)
  num.points <- 1000
  HD <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
              strand = '*',
              hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
    keep.extra.columns = TRUE)
  nlms <- nonlinearFitDist(list(HD), column = 1, verbose = FALSE)

  x=seq(0.1, 10, 0.05)
  y <- predict(nlms[[1]], pred="dens", q = x,
                  dist.name="Weibull2P")
  y1 <- dweibull(x, shape = 0.75, scale = 1)
  # The maximum difference between the "theoretical" and estimated densities
  expect_true(max(abs(round(y, 2) - round(y1, 2))) < 0.1)
})
