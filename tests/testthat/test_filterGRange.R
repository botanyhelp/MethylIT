library(testthat)
library(GenomicRanges)
library(MethylIT)

context("filterGRange tests")

test_that("filterGRange test", {
  gr1 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = c("+","-","+","*","."), mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)
  expected.answer <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:12, end = 11:12,
               strand = c("+","-"), mC = 1, uC = 1:2),
    keep.extra.columns = TRUE)
  answer <- filterGRange(gr1, min.coverage = 1, max.coverage = 4,
                             col.names = c(mC = 1, uC = 2), verbose = FALSE)
  expect_equivalent(answer, expected.answer)
})

test_that("filterGRange empty test", {
  gr1 <- GRanges()
  #TODO: should filterGRange should check the parameters
  expect_error(filterGRange(gr1, min.coverage = 1, max.coverage = 4,
               col.names = c(mC = 1, uC = 2), verbose = FALSE))
})
