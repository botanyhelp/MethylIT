#' @rdname validateClass
#' @title Function to validate S3 classes in MethylIT
#' @param LR An object from class 'pDMP' or 'InfDiv'
#' @export
validateClass <- function(LR, ...) UseMethod("validateClass", LR)

#' @rdname validateClass
#' @name validateClass.default
#' @export
validateClass.default <- function(LR, ...){
  warning(paste("'validateClass' does not know how to handle object of class ",
               class(LR),
               " and can only be used on classes 'pDMP', 'InfDiv'"))

}

#' @rdname validateClass
#' @name validateClass.pDMP
#' @export
validateClass.pDMP <- function(LR, ...) UseMethod("validateClass", LR)
validateClass.pDMP <- function(LR) {
   vn <- c("hdiv", "TV", "wprob")
   if (any(!unlist(lapply(LR, function(GR) class(GR) == "GRanges")))) {
       warning("At least one element from 'LR' is not a 'GRanges' object")
       stop("LR is not a valid 'pDMP' object")
   }
   nams <- unlist(lapply(LR, function(GR) {
       ns <- colnames(mcols(GR))
       sum(is.element(vn, ns))
   }))
   if (any(nams != 3)) {
       warning("At least one element from 'LR' has incorrect column names")
       stop("LR is not a valid 'pDMP' object")
   } else invisible(TRUE)
}

#' @rdname validateClass
#' @name validateClass.InfDiv
#' @export
validateClass.InfDiv <- function(LR, ...) UseMethod("validateClass", LR)
validateClass.InfDiv <- function(LR) {
   vn <- c("hdiv", "TV", "wprob")
   if (any(!unlist(lapply(LR, function(GR) class(GR) == "GRanges")))) {
       warning("At least one element from 'LR' is not a 'GRanges' object")
       stop("LR is not a valid 'InfDiv' object")
   }
   nams <- unlist(lapply(LR, function(GR) {
       ns <- colnames(mcols(GR))
       sum(is.element(vn, ns))
   }))
   if (any(nams != 3)) {
       warning("At least one element from 'LR' has incorrect column names")
       stop("LR is not a valid 'InfDiv' object")
   } else invisible(TRUE)
}


