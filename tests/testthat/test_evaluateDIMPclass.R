library(testthat)
library(GenomicRanges)
library(MethylIT)

context("evaluateDIMPclass tests")

test_that("evaluateDIMPclass dummy test", {
  # set.seed(123) # To set a seed for random number generation
  # GRanges object of the reference with methylation levels in
  # its matacolumn
  num.points = 5000
  Ref = makeGRangesFromDataFrame(
    data.frame(chr = '1',
               start = 1:num.points,
               end = 1:num.points,
               strand = '*',
               p1 = rbeta(num.points, shape1 = 1, shape2 = 1.5)),
    keep.extra.columns = TRUE)

  # List of Granges objects of individuals methylation levels
  Indiv <- GRangesList(
    sample11 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 1:num.points,
                 end = 1:num.points,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 1.5, shape2 = 2)),
      keep.extra.columns = TRUE),
    sample12 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 1:num.points,
                 end = 1:num.points,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 1.6, shape2 = 2.1)),
      keep.extra.columns = TRUE),
    sample21 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 1:num.points,
                 end = 1:num.points,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 40, shape2 = 4)),
      keep.extra.columns = TRUE),
    sample22 = makeGRangesFromDataFrame(
      data.frame(chr = '1',
                 start = 1:num.points,
                 end = 1:num.points,
                 strand = '*',
                 p2 = rbeta(num.points, shape1 = 41, shape2 = 4)),
      keep.extra.columns = TRUE))
  # To estimate Hellinger divergence using only the methylation levels.
  HD <- estimateDivergence(ref = Ref, indiv = Indiv, meth.level = TRUE,
                           columns = 1)
  # To perform the nonlinear regression analysis
  nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)

  # Next, the potential signal can be estimated
  PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 4, alpha = 0.05)

  # The cutpoint estimation used to discriminate the signal from the noise
  cutpoints = estimateCutPoint(PS,
                               control.names = c("sample11", "sample12"),
                               treatment.names = c("sample21", "sample22"),
                               div.col = 4, verbose = TRUE)
  # DIMPs are selected using the cupoints
  DIMPs = selectDIMP(PS, div.col = 4, cutpoint = min(cutpoints$cutpoint))

  # Classification of DIMPs into two clases: DIMPS from control and DIMPs from
  # treatment samples and evaluation of the classifier performance (for more
  # details see ?DIMPclassVal).
  conf.mat <- evaluateDIMPclass(DIMPs,
                                column = c(hdiv = FALSE, TV = FALSE,
                                           wprob = TRUE, pos = FALSE),
                                control.names = c("sample11", "sample12"),
                                treatment.names = c("sample21", "sample22"))

  # The accuracy should be > 0.5
  expect_true(conf.mat$conf.mat$overall[1] > 0.6)
})
