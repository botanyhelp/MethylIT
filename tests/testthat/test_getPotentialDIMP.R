library(testthat)
library(GenomicRanges)
library(MethylIT)

context("getPotentialDIMP tests")

test_that("getPotentialDIMP dummy test", {
  num.points <- 1000
  HD <- GRangesList(
    sample1 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
      keep.extra.columns = TRUE))

  nlms <- nonlinearFitDist(HD, column = 1, verbose = FALSE)

  DIMPs = getPotentialDIMP(LR = HD, nlms = nlms, div.col = 1, alpha = 0.05)
  expect_true(length(HD$sample1) > length(DIMPs$sample1))
})


