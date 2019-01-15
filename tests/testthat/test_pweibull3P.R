library(testthat)
library(GenomicRanges)
library(MethylIT)

context("pweibull3P tests")

test_that("if ", {
  num.samples <- 1000
  shape <- 0.75
  scale <- 1
  x <- rweibull(num.samples, shape = shape, scale = scale)
  wei.model <- weibull3P(x)
  y <- pweibull3P(x,
                  shape = as.numeric(wei.model$Estimate[1]),
                  scale = as.numeric(wei.model$Estimate[2]),
                  mu = as.numeric(wei.model$Estimate[3]) )
  #TODO: How test this?? PASS
  TRUE
})


