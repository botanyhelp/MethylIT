#' @rdname unlistX
#' @aliases unlist.pDMP
#' @aliases unlist.InfDiv
#' @title Flatten Lists extended to 'pDMP' and 'InfDiv' objects
#' @description Given a 'pDMP' or 'InfDiv' objects, unlist simplifies it to 
#'     produce a GRanges object which contains all the GRanges components which
#'     occur in 'pDMP' or 'InfDiv' object. 
#' @param x AN object from the class 'pDMP'.
#' @details This is a method to extend unlist generic function to handle 
#'     the 'pDMP' classes of objects
#' @export
#' @examples 
#' gr1 <-GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
#'           strand = "+", score = 5L, GC = 0.45)
#' gr2 <-
#'   GRanges(seqnames = c("chr1", "chr1"),
#'           ranges = IRanges(c(7,13), width = 3),
#'           strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
#' gr3 <-
#'   GRanges(seqnames = c("chr1", "chr2"),
#'           ranges = IRanges(c(1, 4), c(3, 9)),
#'           strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
#' grl <- list("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
#' # unlist(grl) # It does not work
#' class(grl) <-'InfDiv' # A trick
#' unlist(grl) # It works

unlist <- function(x) UseMethod("unlist", x)

#' @name unlist.pDMP
#' @rdname unlistX
#' @export
unlist.pDMP <- function(x) {
   if (!all(sapply(x, is, "GRanges")))
       stop("all elements in 'x' must be GRanges objects")
   x <- suppressWarnings(do.call("c", unname(x)))
   return(x)
}

#' @name unlist.InfDiv
#' @rdname unlistX
#' @export
unlist.InfDiv <- function(x) {
  if (!all(sapply(x, is, "GRanges")))
    stop("all elements in 'x' must be GRanges objects")
  x <- suppressWarnings(do.call("c", unname(x)))
  return(x)
}


