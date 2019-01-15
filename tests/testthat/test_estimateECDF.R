library(testthat)
library(MethylIT)
context("estimateECDF tests")


test_that("estimateECDF test", {
    x = sample(1:500, 50, replace=TRUE)
    estimateECDF(x, npoints = 15)
    TRUE
})


