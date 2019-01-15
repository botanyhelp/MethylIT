library(testthat)
library(GenomicRanges)
library(MethylIT)

context("uniqueGRanges tests")

test_that("uniqueGRanges dummy test", {
  dfChr1 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
                       strand = c("+","-","+","*","."), score = 1:5)
  dfChr2 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
                       strand = c("+","-","+","*","."), score = 1:5)
  dfChr3 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
                       strand = c("+","-","+","*","."), score = 1:5)

  gr1 <- makeGRangesFromDataFrame(dfChr1, keep.extra.columns = TRUE)
  gr2 <- makeGRangesFromDataFrame(dfChr2, keep.extra.columns = TRUE)
  gr3 <- makeGRangesFromDataFrame(dfChr3, keep.extra.columns = TRUE)

  grList <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)

  r1 <- uniqueGRanges(grList)
  res.expected <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = c("+","-","+","*","."),
               score = 1:5, score.1 = 1:5, score.2 = 1:5),
    keep.extra.columns = TRUE)

  expect_equivalent(r1, res.expected)
    })
