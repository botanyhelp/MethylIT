#' @rdname GeneUpDownStream
#' @title Get Genes plus Up and Down Stream Regions
#' @description Given a genes region or genomic region (GR), this function
#'     yields the GR plus the especified amount of DNA bases upstream and
#'     downstream the GR.
#' @details Users can select whether to request only upstream, only downstream,
#'    or both, upstream and downstream.
#' @param GR A \code{\link[GenomicRanges]{GRanges-class}} object containing the
#'     ranges of the genes or genomic regions to be extended upstream/downstream
#' @param upstream Integer (Default: 0). The amount of DNA bases (bps) upstream
#'     of the GR.
#' @param downstream Integer (Default: 0). The amount of DNA bases (bps)
#'     downstream of the GR.
#' @param onlyUP Logic (Default: FALSE). If TRUE returns the region upstream the
#'     GR.
#' @param onlyDown Logic (Default: FALSE). If TRUE returns the region downstream
#'     the GR.
#' @examples
#' starts = c(65419, 450703, 923928, 944204)
#' ends = c(71585, 451697, 944581, 959309)
#' chrs = c(rep("chr1", 2), rep("chr2", 2))
#' gr = makeGRangesFromDataFrame(
#'         data.frame(seqnames = chrs, start = starts, end = ends,
#'                 strand = c("+", "-", "+", "-"),
#'                 genes = c("A", "B", "C", "D")), keep.extra.columns = TRUE)
#'
#' gr1 = GeneUpDownStream(GR = gr, upstream = 2000, downstream = 1000)

#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @export

GeneUpDownStream <- function(GR, upstream=0, downstream=0, onlyUP=FALSE,
                               onlyDown=FALSE)
{
   if (upstream > 0 && !onlyUP && !onlyDown) {
       strands <- as.character( strand( GR ) )
       starts <- start( GR )
       ends <- end( GR )
       chrs <- seqnames( GR )

       ind <- which( strands == "+" )
       starts[ ind ] <- starts[ ind ] - upstream

       ind <- which( strands == "-" )
       ends[ ind ] <- ends[ ind ] + upstream

       GR.up <- GRanges( seqnames = chrs,
                       ranges = IRanges( start = starts, end = ends ),
                       strand = strands )
       mcols(GR.up) <- mcols(GR)
       GR <- GR.up; rm(GR.up); gc()
   }

   if (downstream > 0 && !onlyUP && !onlyDown) {
       strands <- as.character( strand(GR) )
       starts <- start(GR)
       ends <- end(GR)
       chrs <- seqnames(GR)

       ind <- which( strands == "+" )
       ends[ ind ] <- ends[ ind ] + downstream

       ind <- which( strands == "-" )
       starts[ ind ] <- starts[ ind ] - downstream

       GR.down <- GRanges( seqnames = chrs,
                       ranges = IRanges( start = starts, end = ends ),
                       strand = strands )
       mcols(GR.down) <- mcols(GR)
       GR <- GR.down; rm(GR.down); gc()
   }

   if (onlyUP) {
       if (onlyDown) stop("* If onlyUP is TRUE, then onlyDown must be FALSE")
       strands <- as.character( strand( GR ) )
       starts <- start( GR )
       ends <- end( GR )
       chrs <- seqnames( GR )
       ind <- which( strands == "+" )

       ends[ ind ] <- starts[ ind ] - 1
       starts[ ind ] <- starts[ ind ] - upstream

       ind <- which( strands == "-" )
       starts[ ind ] <- ends[ ind ] + 1
       ends[ ind ] <- ends[ ind ] + upstream

       GR.up <- GRanges( seqnames = chrs,
                        ranges = IRanges( start = starts, end = ends ),
                        strand = strands )
       mcols(GR.up) <- mcols(GR)
       GR <- GR.up; rm(GR.up); gc()
   }

   if (onlyDown) {
       if (onlyUP) stop("* If onlyDown is TRUE, then onlyUP must be FALSE")
       strands <- as.character( strand(GR) )
       starts <- start(GR)
       ends <- end(GR)
       chrs <- seqnames(GR)

       ind <- which( strands == "+" )
       starts[ ind ] <- ends[ ind ] + 1
       ends[ ind ] <- ends[ ind ] + downstream

       ind <- which( strands == "-" )
       ends[ ind ] <- starts[ ind ] - 1
       starts[ ind ] <- starts[ ind ] - downstream

       GR.down <- GRanges(seqnames = chrs,
                       ranges = IRanges( start = starts, end = ends ),
                       strand = strands )
       mcols(GR.down) <- mcols(GR)
       GR <- GR.down; rm(GR.down); gc()
   }
   return(GR)
}
