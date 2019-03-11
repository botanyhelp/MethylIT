#' @rdname predict.LogisticR
#' @name predict.LogisticR
#' @aliases predict.LogisticR
#' @title Predict function for 'LogisticR' method
#' @description Predict using a PCA-LDA model built with function 'LogisticR'
#' @details This function is specific for predictions based on a logistic model
#'     given by function 'evaluateDIMPclass'. A logistic model obtained with
#'     'glm' regression can be used directly with function 'predict' from
#'     'stats' package.
#' @param object To use with function 'predict'. A 'glm' object from a logistic
#'     model containing a list of two objects: 1) an object of class inheriting
#'     from "lda" and 2) an object of class inheriting from "prcomp".
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction
#' @param type To use with function 'predict'. . The type of prediction
#'     required. The default is "all" given by function 'predict.lda' from MASS
#'     package: 'class', 'posterior', and 'scores' (see ?predict.lda).
#' @param ... Not in use.
#'
#' @return A character vector of prediction classes or a numeric vector of
#'   probabilities or a list containing the two vectors: prediction classes
#'   and probabilities.
#'
#' @export
predict.LogisticR <- function(object, ...) UseMethod("predict")
predict.LogisticR <- function(object, newdata,
                              type=c("class", "probabilities", "all"), ...) {
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
   pred <- predict.glm(object, newdata=newdata[, vn], type="response")
   PredClass <- rep( "CT", length(pred))
   PredClass[pred > 0.5] <- "TT"
   pred <- switch(type[1], class=PredClass, probabilities=pred,
               all=list(class=PredClass, probabilities=pred))
   return(pred)
}
