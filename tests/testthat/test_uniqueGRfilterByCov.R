library(testthat)
library(GenomicRanges)
library(MethylIT)

context("uniqueGRfilterByCov tests")

test_that("uniqueGRfilterByCov with overlapped rows", {
   df1 <- data.frame(chr = "chr1", start = 11:16, end = 11:16,
                     mC = c(2,10,7,9,1,10), uC = c(30,20,4,8,0,10))
   df2 <- data.frame(chr = "chr1", start = 12:18, end = 12:18,
                     mC2 = 1:7, uC2 = 0:6)
   gr1 <- makeGRangesFromDataFrame(df1, keep.extra.columns = TRUE)
   gr2 <- makeGRangesFromDataFrame(df2, keep.extra.columns = TRUE)
   r1 <- uniqueGRfilterByCov(gr1, gr2, min.meth = 1, ignore.strand = TRUE)
   expect_true(length(r1) == 2)
})

