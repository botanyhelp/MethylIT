library(testthat)
library(GenomicRanges)
library(MethylIT)

context("pggamma tests")

test_that("pggamma test", {
  q <- (1:9)/10
  p.stats <- pgamma(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
  p.MethylIT <- pggamma(q, alpha = 1, scale = 1, mu = 0,
                        psi = 1, lower.tail = TRUE, log.p = FALSE)
  expect_equivalent(p.MethylIT, p.stats)
})
