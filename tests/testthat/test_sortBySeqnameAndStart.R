library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT tests")

test_that("sortBySeqnameAndStart of empty GRange is an empty GRange", {
  gr1 <- GRanges()
  gr2 <- sortBySeqnameAndStart(gr1)
  expect_equal(gr1, gr2)
})

test_that("sortBySeqnameAndStart same chr and start in different start order", {
  x1 <- c("chr1:1-1", "chr1:2-2")
  x2 <- c("chr1:2-2", "chr1:1-1")
  gr1 = sortBySeqnameAndStart(as(x1, "GRanges"))
  gr2 = sortBySeqnameAndStart(as(x2, "GRanges"))
  expect_identical(gr1, gr2)
})

test_that("sortBySeqnameAndStart same chr and start in different chr order", {
  x1 <- c("chr1:1-1", "chr2:1-1")
  x2 <- c("chr2:1-1", "chr1:1-1")
  gr1 = as(x1, "GRanges")
  gr2 = sortBySeqnameAndStart(as(x2, "GRanges"))
  expect_identical(gr1, gr2)
})

test_that("sortBySeqnameAndStart given chr1:1,2,3 to chr1,2,3 ", {
  # TODO: test all the posible combinations
  # l <- gtools::permutations(n = length(chrs), r = length(chrs), v = chrs)
  TRUE
})


