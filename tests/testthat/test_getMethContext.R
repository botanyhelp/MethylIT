library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT getMethContext tests")

test_that("getMethContext function test", {

    dna<-Biostrings::DNAString(x="CCCTAACGACCCTAACGCTACCCTAAACCTCTGAAT",
     start=1, nchar=NA)
    gr <- getMethContext(chr.seq = dna, chromosome = "1", verbose = TRUE)
    expect_is(gr, "GRanges")
    expect_is(gr$trinucleotide, "factor")
    expect_is(gr$context, "factor")
    expect_is(gr$subcontext, "factor")
    TRUE
    })
