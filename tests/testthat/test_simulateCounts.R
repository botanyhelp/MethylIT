library(testthat)
library(MethylIT)
context("MethylIT simulateCounts tests")

test_that("simulateCounts function test", {
  # *** Simulate samples with expected average of difference of methylation level
  # equal to 0.0427.
  # === The number of cytosine sitesto generate ===
  sites = 100
  # == Set a seed for pseudo-random number generation ===
  set.seed(123)

  # === Simulate samples ===
  ref = simulateCounts(num.samples = 1, sites = sites, alpha = 0.007,
                       beta = 0.5, size = 50, theta = 4.5, sample.ids = "C1")
  treat = simulateCounts(num.samples = 2, sites = sites, alpha = 0.03,
                         beta = 0.5, size = 50, theta = 4.5,
                         sample.ids = c("T1", "T2"))

  #  === Estime Divergences ===
  HD = estimateDivergence(ref = ref$C1, indiv =  treat, Bayesian = TRUE,
                          num.cores = 1L, percentile = 1)
  expect_true(mean(HD$T1$TV) < 0.5)
})
