library(testthat)
library(GenomicRanges)
library(MethylIT)

context("FisherTest tests")

test_that("FisherTest tests", {
   data(HD)
   ## Only the first four cytosine sites from each sample are tested
   hd <- lapply(HD, function(hd) hd[1:4])

   x <- FisherTest(LR = hd, pooling.stat = "sum",
           treatment.names = c("T1", "T2"), tv.cut = NULL,
           pAdjustMethod="BH", pvalCutOff = 0.05, num.cores = 1L,
           verbose=FALSE)
   expect_true(length(x$T3) == 1)
})
