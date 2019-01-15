library(testthat)
library(GenomicRanges)
library(MethylIT)

context("readCounts2GRangesList tests")

test_that("readCounts2GRangesList from a file with 'gz' in the name", {
  # Create a cov file with it's file name including "gz" (tarball extension)
  filename <- "./filegz.cov"
  gr1 = sortBySeqnameAndStart(as(c("chr1:1-1", "chr1:2-2"), "GRanges"))
  write.table(as.data.frame(gr1),
              file = filename, col.names = FALSE, row.names = FALSE,
              quote = FALSE)

  # Load and delete the file
  gr2 <- readCounts2GRangesList(filenames = filename,
                                sample.id = 'test',
                                columns = c(seqnames = 1, start = 2, end = 3),
                                verbose = TRUE)
  file.remove(filename)

  # Test if gr1 == gr2
  expect_true( all(gr1 == gr2$test) )
})

test_that("readCounts2GRangesList with wrong columns argument", {
  filename <- "./file.cov"
  gr1 = sortBySeqnameAndStart(as(c("chr1:1-1", "chr1:2-2"), "GRanges"))
  write.table(as.data.frame(gr1),
              file = filename, col.names = FALSE, row.names = FALSE,
              quote = FALSE)

  # We expect an error!!
  expect_error(
    readCounts2GRangesList(files_names = filename,
                           sample.id = 'test',
                           columns = c(seqnames = 1, start = 2, end = 2),
                           verbose = TRUE))

  # Clean up
  file.remove(filename)
})
