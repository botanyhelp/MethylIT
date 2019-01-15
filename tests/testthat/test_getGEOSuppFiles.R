library(testthat)
library(MethylIT)
context("MethylIT getGEOSuppFiles tests")

test_that("getGEOSuppFiles function test", {
    filenames = getGEOSuppFiles(GEO = "GSM2041700",pattern = "CGmethratio.tab.gz")
    LR13 = readCounts2GRangesList(filenames = filenames,
                            sample.id = c("unknown"),
                            columns = c(seqnames = 1, start = 2,
                                        strand = 3, mC = 4, uC = 5), 
                            remove = TRUE, pattern = "^chr10",
                            chromosome.names = "chr10", verbose = FALSE)
    expect_is(LR13$unknown$mC, "integer")
    file.remove(filenames)
    TRUE
    })
