library(testthat)
library(GenomicRanges)
library(MethylIT)

context("poolFromGRlist tests")

test_that("poolFromGRlist test", {
  gr1 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = '*', mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)
  gr2 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = '*', mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)

  answer <- poolFromGRlist(list(gr1, gr2), stat = 'sum', verbose = FALSE)
  expected.answer <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = '*', mC = 2, uC = seq(2,10,2)),
    keep.extra.columns = TRUE)

  expect_equivalent(answer, expected.answer)
})
