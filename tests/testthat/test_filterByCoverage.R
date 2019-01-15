library(testthat)
library(GenomicRanges)
library(MethylIT)

context("filterByCoverage tests")

test_that("filterByCoverage test", {
  gr1 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = c("+","-","+","*","."), mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)
  expected.answer <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:12, end = 11:12,
               strand = c("+","-"), mC = 1, uC = 1:2),
    keep.extra.columns = TRUE)
  answer <- filterByCoverage(gr1, min.coverage = 1, max.coverage = 4,
                             col.names = c(mC = 1, uC = 2), verbose = FALSE)
  expect_equivalent(answer, expected.answer)
})

test_that("filterByCoverage test GRanges list", {
  gr1 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = c("+","-","+","*","."), mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)
  gr2 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = c("+","-","+","*","."), mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)
  expected.answer <- list(
    makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 11:12, end = 11:12,
                 strand = c("+","-"), mC = 1, uC = 1:2),
      keep.extra.columns = TRUE),
    makeGRangesFromDataFrame(
      data.frame(chr = "chr1", start = 11:12, end = 11:12,
                 strand = c("+","-"), mC = 1, uC = 1:2),
      keep.extra.columns = TRUE))

  answer <- filterByCoverage(list(gr1, gr2), min.coverage = 1, max.coverage = 4,
                             col.names = c(mC = 1, uC = 2), verbose = FALSE)
  expect_equivalent(answer, expected.answer)
})



test_that("filterByCoverage empty GRange test", {
  gr1 <- GRanges()
  #TODO: should filterByCoverage should check the parameters
  #answer <- filterByCoverage(gr1, min.coverage = 1, max.coverage = 4,
  #                           col.names = c(mC = 1, uC = 2), verbose = FALSE)
  #expect_equivalent(answer, expected.answer)
  TRUE
})

