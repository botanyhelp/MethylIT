library(testthat)
library(GenomicRanges)
library(MethylIT)

context("betaBinPosteriors tests")

test_that("function dummy test", {
  res <- MethylIT:::.betaBinPosteriors(1, 4, 0, 0)
  expect_equal(res, 0.25)
})
