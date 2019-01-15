library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT AICmodel tests")

test_that("AICmodel function test", {
    #control_v_treatment_data <- factor(c(rep(0,180),rep(1,20),rep(0,40),rep(1,10)))
    #genotype <-factor(c(rep("AA/Aa",200),rep("aa",50)),levels=c("AA/Aa","aa"))
    #df_data_geno <- data.frame(control_v_treatment_data, genotype)
    #fit <- glm(control_v_treatment_data~genotype, family="binomial", data=df_data_geno) 
    set.seed(77)
    x = runif(100, 1, 5)
    y = 2 * exp(-0.5 * x) + runif(100, 0, 0.1)
    #plot(x, y)
    nlm <- nls(Y ~ a * exp( b * X), data = data.frame(X = x, Y = y),
                start = list( a = 1.5, b = -0.7),
                control = nls.control(maxiter = 10^4, tol = 1e-05),
                algorithm = "port")
    # The estimations of Akaike information criteria given by AIC' function from
    # 'stats' R package and from 'AICmodel' function are equals.
    #AICmodel(nlm) == AIC(nlm)
    #all.equal(AICmodel(nlm), AIC(nlm))
    expect_equal(AICmodel(nlm), AIC(nlm))
    # Now, using residuals from the fitted model:
    res = y - coef(nlm)[1] * exp(coef(nlm)[2] * x)
    expect_equal(AICmodel(residuals = res, np = 2), AIC(nlm))

    answer <- AICmodel(nlm)
    expect_is(answer, "numeric")
    TRUE
    })
