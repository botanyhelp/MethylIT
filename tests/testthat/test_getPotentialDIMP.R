library(testthat)
library(GenomicRanges)
library(MethylIT)

context("getPotentialDIMP tests")

test_that("getPotentialDIMP dummy test", {
   data(HD, gof)
   PS <- getPotentialDIMP(LR = HD, nlms = gof$nlms, dist.name = gof$bestModel,
                           div.col = 9L, alpha = 0.05)
   expect_true(length(HD$T1) > length(PS$T1))
})


