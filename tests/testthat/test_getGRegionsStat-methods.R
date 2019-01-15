library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT getGRegionsStat-methods tests")

test_that("getGRegionsStat-methods.R function test", {
    gr = GRanges(seqnames = Rle( c("chr1", "chr2", "chr3", "chr4"),
             c(5, 5, 5, 5)),
             ranges = IRanges(start = 1:20, end = 1:20),
             strand = rep(c("+", "-"), 10),
             GC = seq(1, 0, length = 20))
    grs = getGRegionsStat(gr, win.size = 4, step.size = 4)
    expect_is(grs, "GRanges")
    # Selecting the positive strand
    grs = getGRegionsStat(gr, win.size = 4, step.size = 4, select.strand = "+")
    expect_is(grs$statistic, "numeric")
    # Selecting the negative strand
    grs = getGRegionsStat(gr, win.size = 4, step.size = 4, select.strand = "-")
    expect_is(grs$statistic, "numeric")
    TRUE
    })
