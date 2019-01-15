library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT predictDIMPclass tests")

test_that("predictDIMPclass function test", {
  num.points <- 5000
  set.seed(123)
  hdiv11 = rweibull(1:num.points, shape = 0.45, scale = 1.2)
  wprob11 = pweibull(hdiv11,
                     shape = 0.45,
                     scale = 1.2,
                     lower.tail = FALSE)
  hdiv12 = rweibull(1:num.points, shape = 0.45, scale = 1.2)
  wprob12 = pweibull(hdiv12,
                     shape = 0.45,
                     scale = 1.2,
                     lower.tail = FALSE)
  hdiv21 = rweibull(1:num.points, shape = 0.6, scale = 1.02)
  wprob21 = pweibull(hdiv21,
                     shape = 0.6,
                     scale = 1.02,
                     lower.tail = FALSE)
  hdiv22 = rweibull(1:num.points, shape = 0.61, scale = 1.02)
  wprob22 = pweibull(hdiv22,
                     shape = 0.61,
                     scale = 1.02,
                     lower.tail = FALSE)
  #' Potential signal
  PS <- GRangesList(
    sample11 = makeGRangesFromDataFrame(
      data.frame(
        chr = "chr1",
        start = 1:num.points,
        end = 1:num.points,
        strand = '*',
        hdiv = hdiv11,
        wprob = wprob11
      ),
      keep.extra.columns = TRUE
    ),
    sample12 = makeGRangesFromDataFrame(
      data.frame(
        chr = "chr1",
        start = 1:num.points,
        end = 1:num.points,
        strand = '*',
        hdiv = hdiv12,
        wprob = wprob12
      ),
      keep.extra.columns = TRUE
    ),
    sample21 = makeGRangesFromDataFrame(
      data.frame(
        chr = "chr1",
        start = 1:num.points,
        end = 1:num.points,
        strand = '*',
        hdiv = hdiv21,
        wprob = wprob21
      ),
      keep.extra.columns = TRUE
    ),
    sample22 = makeGRangesFromDataFrame(
      data.frame(
        chr = "chr1",
        start = 1:num.points,
        end = 1:num.points,
        strand = '*',
        hdiv = hdiv22,
        wprob = wprob22
      ),
      keep.extra.columns = TRUE
    )
  )
  cutpoint = 5.76
  DIMPs = selectDIMP(PS, div.col = 1, cutpoint = cutpoint)

  #' A classification model can be fitted as follow:
  conf.mat <- evaluateDIMPclass(
    DIMPs,
    column = c(
      hdiv = TRUE,
      TV = FALSE,
      wprob = FALSE,
      pos = FALSE
    ),
    interaction = "wprob:hdiv",
    control.names = c("sample11", "sample12"),
    treatment.names = c("sample21", "sample22")
  )
  # Now predictions of DIMP for control and treament can be obtained
  pred = predictDIMPclass(
    LR = DIMPs,
    model = conf.mat$model,
    conf.matrix = TRUE,
    control.names = c("sample11", "sample12"),
    treatment.names = c("sample21", "sample22")
  )

  expect_is(pred, "list")
  expect_is(pred$accuracy, "numeric")
  expect_is(pred$conf.mat, "table")
})
