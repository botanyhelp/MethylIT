library(testthat)
library(GenomicRanges)
library(MethylIT)

context("estimateDivergence tests")

test_that("estimateDivergence test", {
  num.samples = 250
  x <- data.frame(chr = "chr1", start = 1:num.samples, end = 1:num.samples,
                  strand = '*',
                  mC = rnbinom(size = num.samples, mu = 4, n = 500),
                  uC = rnbinom(size = num.samples, mu = 4, n = 500))
  x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  HD <- estimateDivergence(ref = x, indiv = list(x))
  expect_true(all(data.frame(HD)$hdiv == 0))
})
