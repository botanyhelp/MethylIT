library(testthat)
library(GenomicRanges)
library(MethylIT)

context("ksTest tests")

test_that("ksTest test with good estimated parameters", {
  num.samples <- 1000
  x <- rweibull(num.samples, shape = 1.001, scale = 1.001)
  ans <- ksTest(x, pars = c(shape = 1, scale = 1), sample.size = 100)
  # very relaxed test!
  expect_true(ans$p.value > 0.90)
})

test_that("ksTest test with bad estimated parameters", {
  num.samples <- 1000
  x <- rweibull(num.samples, shape = 1.001, scale = 1.001)
  ans <- ksTest(x, pars = c(shape = 0.5, scale = 0.5), sample.size = 100)
  expect_false(ans$p.value > 0.95)
})
