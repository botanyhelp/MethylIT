#' @rdname Extract
#' @title Subscript operator for objects from "InfDiv" and "pDMP" class
#' @param x A list-like or vector-like object
#' @param i Index specifying elements to extract or replace. Index can be a
#'     numeric or character vector. See ?base::\code{\link[base]{Extract}} for
#'     a description of these arguments.
#' @description \eqn{x[i]} returns an object of the length as \eqn{i} preserving
#'     the class from the original object \eqn{x}.
#' @return Same as in ?base::\code{\link[base]{Extract}}, but preserving
#'     the class from the original object \eqn{x}.
#' @seealso base::\code{\link[base]{Extract}}
#' @importFrom S4Vectors extractROWS

#' @name "[.InfDiv"
#' @rdname Extract
#' @export
"[.InfDiv" <- function(x, i) {
  cl <- class(x)
  x <- .subset(x, i)
  attr(x, "class") <- cl
  return(x)
}

#' @name "[.pDMP"
#' @rdname Extract
#' @export
"[.pDMP" <- function(x, i) {
  cl <- class(x)
  x <- .subset(x, i)
  attr(x, "class") <- cl
  return(x)
}

