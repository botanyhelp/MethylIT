library(testthat)
library(GenomicRanges)
library(MethylIT)

context("pcaLDA tests")

test_that("pcaLDA dummy test", {
  data(iris)
  ld = pcaLDA(formula = Species ~., data = iris, n.pc = 1, max.pc = 2,
              scale = TRUE, center = TRUE)
  set.seed(123)
  idx = sample.int(150, 40)
  newdata = iris[idx, 1:4]
  newdata.prediction = predict(ld, newdata = newdata)
  # The confusion matrix
  x = data.frame(TRUE.class = iris$Species[idx],
                 PRED.class = newdata.prediction$class)
  # The accuracy should be > 0.5
  expect_true(all(diag(table(x)) == c(10,13,13)))
})
