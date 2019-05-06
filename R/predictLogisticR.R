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
#'     prediction. Optionally, an object from class "GRanges", a list of
#'     GRanges, "pDMP" or "InfDIv", in which to look for variables with which to
#'     predict. If omitted, the fitted linear predictors are used.
#' @param type To use with function 'predict'. . The type of prediction
#'     required. The default is "class". Possible outputs are:
#'     'class', 'posterior', and 'scores' (see ?predict.lda).
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param ... Not in use.
#'
#' @return A character vector of prediction classes or a numeric vector of
#'     probabilities or a list containing the original 'newdata' with two
#'     columns added in the meta-columns: prediction classes and 'posterior'
#'     probabilities.
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @export
predict.LogisticR <- function(object, ...) UseMethod("predict")
predict.LogisticR <- function(object, newdata = NULL,
                              type=c("class", "posterior", "all"),
                              num.cores = 1, tasks = 0L, ...) {
   if (!is.element(type[1], c("class", "posterior", "all")))
       stop("The type setting '", type[1], "' does not exist")
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
       v <- v[na.omit(match(vn, v))]
       inter <- unlist(lapply(grep("[:]", vn),
                           function(k) strsplit(vn[k], split = ":")[[1]]))
       vn <- union(v, inter)

       if (is.list(newdata)) {
           dt <- lapply(newdata, function(x) {
                           df <- data.frame(hdiv=x$hdiv, TV=x$TV,
                                           logP=log10(x$wprob + 2.2e-308))
                           if (is.element("pos", vn)) df$pos = position(x)
                           df <- df[, vn]
                           return(df)
           })
       }
       else {
           dt <- cbind(hdiv=newdata$hdiv, TV=newdata$TV,
                       logP=log10(newdata$wprob + 2.2e-308))
           if (is.element("pos", vn)) dt$pos = position(newdata)
           dt <- dt[, vn]
       }
   }
   # ---------------------------------------------------------------------#
   object <- structure(object, class=c("glm", "lm"))
   if (is.null(newdata)) dt <- NULL

   if (num.cores > 1 && is.list(dt)) {
       if (Sys.info()['sysname'] == "Linux") {
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")
       newdata <- bplapply(1:length(dt),
                       function(k) {
                           p <- predict.glm(object, newdata = d[[k]],
                                           type="response")
                           PredClass <- rep( "CT", length(p))
                           PredClass[p > 0.5] <- "TT"
                           newdata$PredClass <- PredClass
                           newdata$posterior <- p
                           return(newdata)
                       },
               BPPARAM=bpparam)
   }
   else {
       # === To not depend on BiocParallel package
       if (num.cores == 1 && is.list(dt)) {
           newdata <- lapply(1:length(dt),
                               function(k) {
                                   p <- predict.glm(object, newdata = d[[k]],
                                                   type="response")
                                   PredClass <- rep( "CT", length(p))
                                   PredClass[p > 0.5] <- "TT"
                                   newdata$PredClass <- PredClass
                                   newdata$posterior <- p
                                   return(newdata)
                               }
                       )
       }
   }
   if (!is.list(dt)) {
       pred <- predict.glm(object, newdata = dt, type="response")
       newdata$PredClass <- rep( "CT", length(pred))
       newdata$PredClass[pred > 0.5] <- "TT"
       newdata$posterior <- pred
   }
   pred <- switch(type[1],
                   class = newdata$PredClass,
                   posterior = newdata$pred,
                   all = newdata)
   return(pred)
}
