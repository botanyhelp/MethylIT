library(testthat)
library(MethylIT)
context("MethylIT classPerform tests")

test_that("classPerform function test", {
  # load simulated data of potential methylated signal
  data(sim_ps)

  cp = classPerform(LR = PS, min.tv = 0.25, tv.cut = 0.4,
                    cutoff = 68.7, tv.col = 7L, div.col = 9, stat = 0)
  expect_true(cp$FDR < 0.1)
})
