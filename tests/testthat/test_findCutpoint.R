library(testthat)
library(MethylIT)
context("MethylIT findCutpoint tests")

test_that("findCutpoint function test", {
  # load simulated data of potential methylated signal
  data(sim_ps)

  # Vector of cutoff values
  cuts = c(2, 5, 10, 15, 18, 20, 21, 22, 25, 27, 30, 35, 40,
           45, 50, 55)
  # # === To find the cutpoint that maximize the accuracy ===
  pre.cut.acc = findCutpoint(LR = PS, min.tv = 0.25, tv.cut = 0.5,
                             predcuts = cuts, tv.col = 7L, div.col = 9,
                             stat = 1, num.cores = 15)
  expect_true(pre.cut.acc$FDR < 0.5)
})
