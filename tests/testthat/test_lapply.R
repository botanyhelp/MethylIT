library(testthat)
library(MethylIT)
context("MethylIT lapply tests")

test_that("lapply function test", {
  x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
  class(x) <- "nice"
  expect_true(class(lapply(x, mean, keep.attr = TRUE)) == "nice")
})
