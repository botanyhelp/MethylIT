#' @rdname evaluateModel
#'
#' @title Evaluate a model using the Akaike information criterion (AIC)
#' @description Evaluate a glm object using the Akaike information criterion
#'     (AIC)
#'
#' @param model glm object
#' @param test A character string matching one of "Wald" or "LRT". If test =
#'   "Wald", then the p-value of the Wald test for the coefficient of the
#'   independent variable (\emph{treatment group}) will be reported p-value. If
#'   test = "LRT", then the p-value from a likelihood ratio test given by
#'   \code{\link[stats]{anova}} function from \emph{stats} packages will be the
#'   reported p-value for the group comparison.
#' @return AIC value
#' @keywords internal
.evaluateModel <- function(model, test = c("Wald", "LRT")) {
   if (!inherits(model, "try-error")) {
       cfs <- coef(summary(model))
       log2FC <- cfs[2, 1]
       if (test[1] == "LRT") {
           anov <- try(suppressWarnings(anova(model, test = "Chisq")),
                       silent=TRUE)
           if (!inherits(anov, "try-error")) coef.pval <- anov$`Pr(>Chi)`[2]
           else coef.pval <- cfs[2, 4]
       } else coef.pval <- cfs[2, 4]
       if (model$family$family == "quasipoisson") {
           AICs <- .AICquasiPoisson(model)
       } else {
           AICs <- AIC(model)
       }
       if (coef.pval < 0.05) {
           res <- list(Eval=TRUE, log2FC=log2FC, coef.pval=coef.pval, AIC=AICs,
                   mdl=model)
       } else {
           res <- list(Eval=FALSE, AIC=1e16)
       }
   } else {
       res <- list(Eval=FALSE, AIC=1e16)
   }
   return(res)
}
