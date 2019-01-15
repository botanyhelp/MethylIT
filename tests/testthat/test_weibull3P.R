library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT tests")

test_that("weibull3P with 2 parameters test", {
  num.samples <- 1000
  shape <- 0.75
  scale <- 1
  x <- rweibull(num.samples, shape = shape, scale = scale)
  res <- weibull3P(x)

  ## TODO: Consider to test the results with Adj.R.Square, rho, and R.Cross.val
  #3 or use two different datasets to cross the results
  ##y <- rweibull(num.samples,
  ##              shape = as.numeric(res$Estimate[1]),
  ##              scale = as.numeric(res$Estimate[2]))
  ##ks.test(x, y, exact = FALSE)
  ##wilcox.test(x, y)

  # less restringent R.Cross.val requirement
  expect_true(as.numeric(res$R.Cross.val[1]) > 0.95)
})
