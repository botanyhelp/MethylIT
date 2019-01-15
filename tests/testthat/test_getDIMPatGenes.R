library(testthat)
library(GenomicRanges)
library(MethylIT)
library(rtracklayer)


context("getDIMPatGenes tests")

test_that("getDIMPatGenes dummy test", {
  genes = GRanges(seqnames = "1",
                  ranges = IRanges(start = c(3631, 6788, 11649),
                                   end = c(5899, 9130, 13714)),
                  strand = c("+", "-", "-"))
  mcols(genes) <- data.frame(gene_id = c("AT1G01010", "AT1G01020",
                                         "AT1G01030"))
  set.seed(135)
  num.points = 5001
  end.point = 8600
  Ref = makeGRangesFromDataFrame(
    data.frame(chr = '1',
               start = 3600:end.point,
               end = 3600:end.point,
               strand = '*',
               p1 = rbeta(num.points, shape1 = 1, shape2 = 1.5)),
    keep.extra.columns = TRUE)

  # List of Granges objects of individuals methylation levels
  Indiv <- GRangesList(
    sample11 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 3600:end.point,
                 end = 3600:end.point,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 1.5, shape2 = 2)),
      keep.extra.columns = TRUE),
    sample12 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 3600:end.point,
                 end = 3600:end.point,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 1.6, shape2 = 2.1)),
      keep.extra.columns = TRUE),
    sample21 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 3600:end.point,
                 end = 3600:end.point,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 33, shape2 = 12)),
      keep.extra.columns = TRUE),
    sample22 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 3600:end.point,
                 end = 3600:end.point,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 35, shape2 = 12)),
      keep.extra.columns = TRUE))
  # To estimate Hellinger divergence using only the methylation levels.
  HD <- estimateDivergence(ref = Ref, indiv = Indiv, meth.level = TRUE,
                           columns = 1, verbose = FALSE)
  # To perform the nonlinear regression analysis
  nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)

  # Next, the potential signal can be estimated
  PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 4, alpha = 0.05)

  # The cutpoint estimation used to discriminate the signal from the noise
  cutpoints = estimateCutPoint(PS,
                               control.names = c("sample11", "sample12"),
                               treatment.names = c("sample21", "sample22"),
                               div.col = 4, verbose = FALSE)
  # DIMPs are selected using the cupoints
  DIMPs = selectDIMP(PS, div.col = 4, cutpoint = min(cutpoints$cutpoint))
  DIMR = getDIMPatGenes(GR = DIMPs$sample22, GENES = genes)
  expect_true(length(DIMR) > 0)
})


