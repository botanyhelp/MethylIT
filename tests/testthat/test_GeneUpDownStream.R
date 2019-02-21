library(testthat)
library(MethylIT)

context("GeneUpDownStream tests")

test_that("GeneUpDownStream test", {
  starts = c(65419, 450703, 923928, 944204)
  ends = c(71585, 451697, 944581, 959309)
  chrs = c(rep("chr1", 2), rep("chr2", 2))
  gr = makeGRangesFromDataFrame(
    data.frame(seqnames = chrs, start = starts, end = ends,
               strand = c("+", "-", "+", "-"),
               genes = c("A", "B", "C", "D")), keep.extra.columns = TRUE)

  gr1 = GeneUpDownStream(GR = gr, upstream = 2000, downstream = 1000)
  expect_true(all((start(gr[c(1,3)]) - 2000) == start(gr1[c(1,3)])))
})
