library(testthat)
library(GenomicRanges)
library(MethylIT)

context("estimateCutPoint tests")

test_that("estimateCutPoint dummy tests", {
  set.seed(123)
  num.points <- 1000

  # A list of GRanges objects with simulated Hellinger divergences in their
  # metacolumns.
  HD <- GRangesList(
    sample11 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 1.3, scale = 1.02)),
      keep.extra.columns = TRUE),
    sample12 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 1.2, scale = 1.02)),
      keep.extra.columns = TRUE),
    sample21 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 1.45, scale = 1.5)),
      keep.extra.columns = TRUE),
    sample22 = makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
                 strand = '*',
                 hdiv = rweibull(1:num.points, shape = 1.75, scale = 1.6)),
      keep.extra.columns = TRUE))

  # Nonlinear fit of Weiblul distribution
  nlms <- nonlinearFitDist(HD, column = 1, verbose = FALSE)

  # Estimation of the potential signal and cutpoints
  PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 1, alpha = 0.05)
  cutpoints = estimateCutPoint(PS, control.names = c("sample11", "sample12"),
                                 treatment.names = c("sample21", "sample22"),
                                 div.col = 1, verbose = FALSE)
  # The AUC should be > 0.5
  expect_true(cutpoints$auc$sample11[1] > 0.5)
})
