#' @rdname predict.LogisticR
#' @name predict.LogisticR
#' @aliases predict.LogisticR
#' @title Predict function for logistic regression model from 'LogisticR' class
#' @description Predict using a logistic model obtained from the output of
#'     function \code{\link{evaluateDIMPclass}}.
#' @details This function is specific for predictions based on a logistic model
#'     given by function \code{\link{evaluateDIMPclass}}. A logistic model is
#'     obtained with 'glm' regression can be used directly with function
#'     'predict' from 'stats' package.
#' @param object To use with function 'predict'. An object from 'LogisticR'
#'     class. A logistic model given by function
#'     \code{\link{evaluateDIMPclass}}.
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction. Optionally, a data frame in which to look for variables with
#'     which to predict. If omitted, the fitted linear predictors are used.
#' @param type To use with function 'predict'. . The type of prediction
#'     required. The default is "class". Possible outputs are:
#'     'class', 'posterior', and 'scores' (see ?predict.lda).
#' @param ... Not in use.
#'
#' @return A character vector of prediction classes or a numeric vector of
#'   probabilities or a list containing the two vectors: prediction classes
#'   and 'posteriors' probabilities.
#'
#' @export
predict.LogisticR <- function(object, ...) UseMethod("predict")
predict.LogisticR <- function(object, newdata = NULL,
                              type=c("class", "posteriors", "all"), ...) {
   if (!is.element(type, c("class", "posteriors", "all")))
       stop("The type setting '", type, "' does not exist")
   if (!inherits(object, "LogisticR")) {
       stop("* 'object' must be logistic a model from class 'LogisticR'")
   }
   v <- c("hdiv", "TV", "logP", "pos")
   vn <- setdiff(names(coef(object)),"(Intercept)")
   v <- v[is.element(vn, v)]
   inter <- unlist(lapply(grep("[:]", vn),
                        function(k) strsplit(vn[k], split = ":")[[1]]))
   vn <- union(v, inter)

   object <- structure(object, class=c("glm", "lm"))
   if (!is.null(newdata)) newdata=newdata[, vn]
   pred <- predict.glm(object, newdata=newdata, type="response")
   PredClass <- rep( "CT", length(pred))
   PredClass[pred > 0.5] <- "TT"
   pred <- switch(type[1], class=PredClass, posteriors=pred,
               all=list(class=PredClass, posteriors=pred))
   return(pred)
}
