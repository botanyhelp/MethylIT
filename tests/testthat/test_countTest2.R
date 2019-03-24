library(testthat)
library(GenomicRanges)
library(MethylIT)

context("countTest2 tests")

test_that("countTest2 test", {
  countData <- matrix(1:40,ncol = 4)
  colnames(countData) <- c("A1","A2","B1","B2")
  colData <- data.frame(condition = factor(c("A","A","B","B")),
                        c("A1","A2","B1","B2"),
                        row.names = 2)
  ds <- glmDataSet(counts = countData, colData = colData)
  y <- countTest2(DS=ds, num.cores = 1L, maxGrpCV = c(0.6, 0.5),
                   verbose = FALSE)
  expect_true(all(y$adj.pval < 0.1))
})
