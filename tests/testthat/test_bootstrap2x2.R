library(testthat)
library(GenomicRanges)
library(MethylIT)

context("bootstrap2x2 tests")

test_that("bootstrap2x2 test", {
  set.seed(123)
  TeaTasting = matrix(c(8, 350, 2, 20), nrow = 2,
  dimnames = list(Guess = c("Milk", "Tea"), Truth = c("Milk", "Tea")))

  # Small num.permut for test's speed sake
  bootstrap2x2( TeaTasting, stat = "all", num.permut = 100 )
  TRUE
})
