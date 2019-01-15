library(testthat)
library(MethylIT)

context("MethylIT colorBar tests")

test_that("colorBar function test", {
    # function color.bar output is to generate an image, as this code does:
    #  We could test with lionel's vdiffr:
    #  install.packages("vdiffr")
    #  but a simple test that it returns NULL is all we do here
    z <- matrix(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5)),nrow=5)
    breaks <- seq(min( z, na.rm = TRUE), max(z, na.rm = TRUE), length.out = 100)
    bar.palette <- colorRampPalette(c(rep("cyan",4), "green",rep("yellow", 2),
                                   rep("red", 3), rep("darkblue", 2),
                                   rep("black",2)), bias = 1.1, space = "rgb")
    img <- MethylIT:::.colorBar(d, col = bar.palette(length( breaks ) - 1), breaks = breaks,
            horiz = FALSE, lwd = 1, cex.lab = 2)
    expect_identical(img,NULL)
    })
