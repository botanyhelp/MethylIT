library(testthat)
library(GenomicRanges)
library(MethylIT)

context("countTest tests")

test_that("countTest test", {
  countData <- matrix(1:40,ncol = 4)
  condition <- factor(c("A","A","B","B"))
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
  answer <- countTest( dds, verbose = FALSE )
  # this data set do not fit any model provided in countTest
  expect_true(all(answer$adj.pval < 0.95))
})
