#' @rdname Extract
#' @name "["
#' @title S3 overloading of subscript operator to preserve class
#' @param x A list-like or vector-like object
#' @param i,... Index specifying elements to extract or replace. Indices are
#'     numeric or character vectors. See ?base::\code{\link[base]{Extract}} for
#'     a description of these arguments.
#' @description \eqn{x[i]} returns an object of the length as \eqn{i} preserving
#'     the class from the original object \eqn{x}.
#' @return Same as in ?base::\code{\link[base]{Extract}}, but preserving
#'     the class from the original object \eqn{x}.
#' @seealso base::\code{\link[base]{Extract}}
#' @examples
#' # Create a list
#' x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
#' class(x) <- "nice"
#'
#' # Preserves attributes
#' x[1:2]
#' @export
#'
"[" <- function(x,i, ...) {
  cl <- class(x)
  x <- base::.subset(x,i)
  attr(x, "class") <- cl
  return(x)
}
