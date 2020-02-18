library(testthat)
library(GenomicRanges)
library(MethylIT)

context("FisherTest tests")

test_that("FisherTest tests", {
   data(HD)
   ### --- To get the read counts
   hd <- lapply(HD, function(hd) {
       hd = hd[1:7,3:4]
       colnames(mcols(hd)) <- c("mC", "uC")
       return(hd)
   })
   x <- FisherTest(LR = hd,
                   pooling.stat = NULL,
                   control.names = "C1",
                   treatment.names = "T1",
                   pAdjustMethod="BH",
                   pvalCutOff = 0.05,
                   num.cores = 1L,
                   verbose=FALSE)

   expect_true(length(x[[1]]) == 3)
})
