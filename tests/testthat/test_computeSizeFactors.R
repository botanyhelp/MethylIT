library(testthat)
library(GenomicRanges)
library(DESeq2)
library(MethylIT)
context("MethylIT tests")

test_that("function dummy test", {
  dds <- DESeq2::makeExampleDESeqDataSet(n = 1000, m = 4)
  MethylIT:::.computeSizeFactors(dds)
  TRUE
})
