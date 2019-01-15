library(testthat)
library(GenomicRanges)
library(MethylIT)

context("selectDIMP tests")

test_that("selectDIMP test", {
  num.points <- 1000
  HD <- GRangesList(
    sample1 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
      keep.extra.columns = TRUE),
    sample2 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 0.75, scale = 1)),
      keep.extra.columns = TRUE))

  nlms <- nonlinearFitDist(HD, column = 1, verbose = FALSE)

  PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 1, alpha = 0.05)
  cutpoints = estimateCutPoint(PS, control.names = "sample1",
                               treatment.names = c("sample2"),
                               div.col = 1, verbose = FALSE)
  DIMPs = selectDIMP(PS, div.col = 1, cutpoint = cutpoints$cutpoint$sample1)
  expect_true(length(DIMPs$sample1) < length(PS$sample1))
})

