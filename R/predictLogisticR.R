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
#'   and 'posterior' probabilities.
#'
#' @export
predict.LogisticR <- function(object, ...) UseMethod("predict")
predict.LogisticR <- function(object, newdata = NULL,
                              type=c("class", "posterior", "all"), ...) {
   if (!is.element(type[1], c("class", "posterior", "all")))
       stop("The type setting '", type, "' does not exist")
   if (!inherits(object, "LogisticR")) {
       stop("* 'object' must be logistic a model from class 'LogisticR'")
   }

   # ---This builds data frames from newdata from class 'pDMP' or 'InfDiv' ----#
   if (!is.null(newdata) &&
           (inherits(newdata, 'pDMP') || inherits(newdata, 'InfDiv'))) {
       if (!validateClass(newdata)) {
           stop("newdata is not a valid object from class '",
               class(newdata),"'" )
       }
       position <- function(gr) {
           chrs <- split(gr, seqnames(gr))
           gr <- lapply(chrs, function(grc) {
               x <- start(grc)
               x.min <- min(x)
               x.max <- max(x)
               delta <-  max(c(x.max - x, 1))
               return((x - x.min) / (delta))})
           return(unlist(gr))
       }
       v <- c("hdiv", "TV", "logP", "pos")
       vn <- setdiff(names(coef(object)),"(Intercept)")
       v <- v[is.element(vn, v)]
       inter <- unlist(lapply(grep("[:]", vn),
                           function(k) strsplit(vn[k], split = ":")[[1]]))
       vn <- union(v, inter)

       dt <- data.frame()
       for (k in 1:length(newdata)) {
           dc <- c()
           x <- newdata[[k]]
           dc <- cbind(hdiv=x$hdiv, TV=x$TV, logP=log10(x$wprob + 2.2e-308),
                       pos = position(x))
           dt <- rbind(dt, data.frame(dc))
       }
       dt <- dt[, vn]
   }
   # ---------------------------------------------------------------------#

   object <- structure(object, class=c("glm", "lm"))
   if (!is.null(newdata)) newdata = dt
   pred <- predict.glm(object, newdata=newdata, type="response")
   PredClass <- rep( "CT", length(pred))
   PredClass[pred > 0.5] <- "TT"
   pred <- switch(type[1], class=PredClass, posterior=pred,
               all=list(class=PredClass, posterior=pred))
   return(pred)
}
