library(testthat)
library(GenomicRanges)
library(MethylIT)

context("estimateHellingerDiv tests")

test_that("estimateHellingerDiv using p", {
  p <- c(0.5, 0.5)
  expect_equal(estimateHellingerDiv(p), 0)
})

test_that("estimateHellingerDiv using p < 1", {
  #TODO: handle errors
  TRUE
})

test_that("estimateHellingerDiv using n", {
   p <- c(0, 2)
   n <- c(100, 200)
   estimateHellingerDiv(p = p, n = n)
   #TODO: handle this error
   TRUE
})
