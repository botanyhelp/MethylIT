#' @name getGRegionsStat-methods
#' @rdname getGRegionsStat-methods
#' @title Statistic of Genomic Regions
#' @description A function to estimate the centrality measures of a specified
#'     variable given in GRanges object (a column from the metacolums of the
#'     GRanges object) after split the GRanges object into intervals.
#' @details This function split a Grange object into intervals genomic regions
#'     (GR) of fixed size (as given in function "tileMethylCounts2" R package
#'     methylKit, with small changes). A summarized statistic (mean, median,
#'     geometric mean or sum) is calculated for the specified variable values
#'     from each region. Notice that if win.size == step.size, then
#'     non-overlapping windows are obtained.
#' @param GR A Grange object with the variable of interest in its metacolumn.
#' @param win.size An integer for the size of the windows/regions size of the
#'     intervals of genomics regions.
#' @param step.size Interval at which the regions/windows must be defined
#' @param grfeatures A GRanges object corresponding to an annotated genomic
#'     feature. For example, gene region, transposable elements, exons,
#'     intergenic region, etc. If provided, then parameters 'win.size' and
#'     step.size are ignored and the statistics are estimated for 'grfeatures'.
#' @param stat Statistic used to estimate the summarized value of the variable
#'     of interest in each interval/window. Posible options are: "mean",
#'     geometric mean ("gmean"), "median", "density", and "sum" (default). Here,
#'     we define "density" as the sum of values from the variable of interest
#'     in the given region devided by the length of the region.
#' @param absolute Optional. Logic (default: FALSE). Whether to use the absolute
#'     values of the variable provided
#' @param select.strand Optional. If provided,"+" or "-", then the summarized
#'     statistic is computed only for the specified DNA chain.
#' @param column Integer number denoting the column where the variable of
#'     interest is located in the metacolumn of the GRanges object or an integer
#'     vector of two elements (only if prob = TRUE).
#' @param prob Logic. If TRUE and the variable of interest has values between
#'     zero and 1, then the summarized statistic is comuputed using Fisher's
#'     transformation. If length(column) == 2, say with colums x1 and x2, then
#'     the variable of interest will be p = x1/(x1 + x2). For example, if x1
#'     and x2 are methylated and unmethylated read counts, respectively, then p
#'     is the methylation level.
#' @param entropy Logic. Whether to compute the entropy when prob == TRUE.
#' @param maxgap,minoverlap,type See ?findOverlaps in the IRanges package for a
#'     description of these arguments.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations.
#' @param scaling integer (default 1). Scaling factor to be used when
#'     stat = "density". For example, if scaling = 1000, then density * scaling
#'     denotes the sum of values in 1000 bp.
#' @param logbase A positive number: the base with respect to which logarithms
#      are computed when parameter 'entropy = TRUE' (default: logbase = 2).
#' @param na.rm Logical value. If TRUE, the NA values will be removed
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @return A GRanges object with the new genomic regions and their corresponding
#'     summarized statistic.
#' @examples
#' gr <- GRanges(seqnames = Rle( c("chr1", "chr2", "chr3", "chr4"),
#'             c(5, 5, 5, 5)),
#'             ranges = IRanges(start = 1:20, end = 1:20),
#'             strand = rep(c("+", "-"), 10),
#'             GC = seq(1, 0, length = 20))
#' grs <- getGRegionsStat(gr, win.size = 4, step.size = 4)
#' grs
#'
#' ## Selecting the positive strand
#' grs <- getGRegionsStat(gr, win.size = 4, step.size = 4, select.strand = "+")
#' grs
#'
#' ## Selecting the negative strand
#' grs <- getGRegionsStat(gr, win.size = 4, step.size = 4, select.strand = "-")
#' grs
#'
#' ## Operating over a list of GRanges objects
#' gr2 <- GRanges(seqnames = Rle( c("chr1", "chr2", "chr3", "chr4"),
#'                             c(5, 5, 5, 5)),
#'                 ranges = IRanges(start = 1:20, end = 1:20),
#'                 strand = rep(c("+", "-"), 10),
#'                 GC = runif(20))
#'
#' grs <- getGRegionsStat(list(gr1 = gr, gr2 = gr2), win.size = 4, step.size = 4)
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom stats median
#' @export
#' @author Robersy Sanchez
#'
#' @aliases getGRegionsStat
#' @rdname getGRegionsStat-methods
setGeneric("getGRegionsStat",
           function(GR, win.size=350, step.size=350, grfeatures=NULL,
                 stat=c("sum", "mean", "gmaean", "median", "density"),
                 absolute=FALSE, select.strand=NULL, column=1L,
                 prob=FALSE, entropy=FALSE, maxgap=-1L,
                 minoverlap=0L, scaling=1000L, logbase = 2,
                 type=c("any", "start", "end", "within", "equal"),
                 ignore.strand=FALSE, na.rm=TRUE,
                 num.cores = 1L, tasks = 0)
             standardGeneric("getGRegionsStat"))

#' @aliases getGRegionsStat
#' @rdname getGRegionsStat-methods
setMethod("getGRegionsStat", signature(GR="GRanges"),
           function(GR, win.size, step.size, grfeatures, stat, absolute,
                   select.strand, column, prob, entropy, maxgap, minoverlap,
                   scaling, logbase, type, ignore.strand, na.rm) {
           ## These NULL quiet: no visible binding for global variable 'x2'
           x1 <- x2 <- ent <- statistic <- NULL
           if (class( GR ) != "GRanges" )
               stop( "object must be a GRanges object!")
           if (!is.null(grfeatures) && !inherits(grfeatures,"GRanges"))
           {
               stop("* 'grfeatures', if provided, must be a GRanges object")
           }

           ## === Some functions to use ===
           statist <- function(x, stat, absolute) {
               if (absolute) x = abs(x)
               x <- switch(stat[1],
                       sum=sum(x, na.rm=na.rm),
                       mean=mean(x, na.rm=na.rm),
                       gmean=exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)),
                       median=median(x, na.rm=na.rm),
                       density=sum(x, na.rm=na.rm))
               }

           ## =============================== ##
           if (!is.null(select.strand)) {
               ## possible values "-", "+", NULL
               if (!is.element(select.strand, unique(strand(GR)))) {
                   stop("* The GRanges object does not have strand named ", "'",
                       select.strand, "'")
               }
               idx <- which(as.character(strand(GR)) == select.strand)
               GR <- GR[idx]
           }

           GR <- GR[, column]

           chrs <- as.character(unique(seqnames(GR)))

           ## === If genomic features are not specified ===
           if (is.null(grfeatures)) {
               if (length(GR) < win.size || length(GR) < step.size) {
                   stop("* 'GR'length is lesser of 'win.size' or 'step.size'")
               }
               all.wins = GRanges()
               for (k in 1:length(chrs)) {
                   ## get max length of chromosome
                   max.length <- max(IRanges::end(GR[seqnames(GR) == chrs[k],]))
                   ## get sliding windows
                   numTiles <- floor((max.length -
                                   (win.size - step.size)) / step.size) + 1
                   ranges <- IRanges(start=(1 + 0:(numTiles - 1) * step.size),
                                   width=rep(win.size, numTiles))
                   temp.wins <- GRanges(seqnames=rep(chrs[k], numTiles),
                                   ranges=ranges)
                   all.wins <- suppressWarnings(c(all.wins, temp.wins))
               }

               ## sites of interest inside of the windows
               Hits <- findOverlaps(GR, all.wins, maxgap=maxgap,
                                   minoverlap=minoverlap,
                                   ignore.strand=ignore.strand,
                                   type=type)
               all.wins <- all.wins[subjectHits(Hits)]
               mcols(all.wins) <- mcols(GR[queryHits(Hits)])
               chr <- seqnames(all.wins)

               ## Variable to mark the limits of each GR
               cluster.id <- data.frame(cluster.id=paste(chr, start(all.wins),
                                   end(all.wins), sep = "_"))
               GR <- all.wins; rm(all.wins); gc()
           } else {
               ## sites of interest inside of the windows
               Hits <- findOverlaps(GR, grfeatures, maxgap=maxgap,
                           minoverlap=minoverlap, ignore.strand=ignore.strand,
                           type=type)
               grfeatures <- grfeatures[subjectHits(Hits)]
               mcols(grfeatures) <- mcols(GR[ queryHits(Hits)])
               chr <- seqnames(grfeatures)
               if (class(names(grfeatures)) == "character") {
                   cluster.id <- data.frame(cluster.id=names(grfeatures))
                   names(grfeatures) <- NULL
               } else {
                   cluster.id <- data.frame(cluster.id=paste(chr,
                                                             start(grfeatures),
                                                             end(grfeatures),
                                                             strand(grfeatures),
                                                             sep="_"))
               }
               GR <- grfeatures; rm(grfeatures); gc()
           }

           mcols(GR) <- DataFrame(cluster.id, mcols(GR))
           GR <- data.table(as.data.frame(GR))
           if (length(column) < 2) {
               colnames(GR) <- c("seqnames", "start", "end", "width",
                               "strand", "cluster.id", "statistic")
           } else {
               colnames(GR) <- c("seqnames", "start", "end", "width",
                               "strand", "cluster.id", "x1", "x2")
           }

           if (prob && length(column) < 2) { ## Apply Fisher transformation
               GR$p <- GR$statistic
               GR$statistic <- atanh(GR$statistic)
           }
           grn <- c("seqnames", "start", "end")
           ## Compute statistic for regions
           if (length(column) < 2) {
               if (!prob && !entropy) {
                   GR <- GR[, list(seqnames=unique(seqnames), start=min(start),
                           end=max(end),
                           statistic=statist(statistic, stat, absolute)),
                           by=cluster.id]
                   GR <- data.frame(GR)[, c(grn, "statistic")]
               }
               if (prob && !entropy) {
                   GR <- GR[, list(seqnames=unique(seqnames), start=min(start),
                           end=max(end),
                           stat.prob=statist(statistic, stat, absolute)),
                           by=cluster.id]
                   GR <- data.frame(GR)[, c(grn, "stat.prob")]
               }
               if (prob && entropy) {
                   GR$ent <- shannonEntr(GR$p, logbase=logbase)
                   GR <- GR[ ,list(seqnames=unique(seqnames), start=min(start),
                           end=max(end),
                           stat.prob=statist(statistic, stat, absolute),
                           stat.ent=statist(ent, stat, absolute)),
                           by=cluster.id]
                   GR <- data.frame(GR)[, c(grn, "stat.prob", "stat.ent")]
               }
               if (entropy && !prob) {
                   GR$ent <- shannonEntr(GR$p, logbase=logbase)
                   GR <- GR[, list(seqnames=unique(seqnames), start=min(start),
                           end=max(end),
                           stat.ent=statist(ent, stat, absolute)),
                           by=cluster.id]
                   GR <- data.frame(GR)[, c(grn, "stat.ent")]
               }
           } else {
               if (prob || entropy) {
                   GR$p <- GR$x1/(GR$x1 + GR$x2)
                   GR$p[ is.na(GR$p) ] <- 0
                   if (prob && !entropy) {
                       GR <- GR[, list(seqnames=unique(seqnames),
                               start=min(start),
                               end=max(end),
                               x1=statist(x1, stat, absolute),
                               x2=statist(x2, stat, absolute)),
                               by=cluster.id]
                       GR$stat.prob <- GR$x1 / (GR$x1 + GR$x2)
                       GR <- data.frame(GR)[, c(grn, "stat.prob")]
                   }
                   if (prob && entropy) {
                       GR$ent <- shannonEntr(GR$p, logbase=logbase)
                       GR <- GR[, list(seqnames=unique(seqnames),
                           start=min(start),
                           end=max(end),
                           x1=statist(x1, stat, absolute),
                           x2=statist(x2, stat, absolute),
                           stat.ent=statist(ent, stat, absolute)),
                           by=cluster.id]
                       GR$stat.prob <- GR$x1 / (GR$x1 + GR$x2)
                       GR <- data.frame(GR)[, c(grn, "stat.prob", "stat.ent")]
                   }
                   if (entropy && !prob) {
                       GR$ent <- shannonEntr(GR$p, logbase=logbase)
                       GR <- GR[, list(seqnames=unique(seqnames),
                           start=min(start),
                           end=max(end),
                           stat.ent=statist(ent, stat, absolute)),
                           by=cluster.id]
                       GR <- data.frame(GR)[, c(grn, "stat.ent")]
                   }
               }
           }

           if (prob && length(column) < 2) {
               ## Apply inverse of Fisher transformation
               GR$stat.prob <- tanh(GR$stat.prob)
           }

           if (!is.null(select.strand)) GR$strand <- select.strand
           GR <- makeGRangesFromDataFrame(GR, keep.extra.columns=TRUE)
           if (stat == "density" && !prob && !entropy) {
             widths=width(GR)
             GR$statistic <- (scaling * GR$statistic/widths)
           }
           return(GR)
          }
)

#' @aliases getGRegionsStat, list-method
#' @rdname getGRegionsStat-methods
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#'
setMethod("getGRegionsStat", signature(GR="list"),
           function(GR, win.size=350, step.size=350, grfeatures=NULL,
                   stat=c("sum", "mean", "gmaean", "median", "density"),
                   absolute=FALSE, select.strand=NULL, column=1L,
                   prob=FALSE, entropy=FALSE, maxgap=-1L, minoverlap=0L,
                   scaling=1000L, logbase = 2,
                   type=c("any", "start", "end", "within", "equal"),
                   ignore.strand=FALSE, na.rm=TRUE, num.cores = 1L, tasks = 0) {
           if (inherits(GR, "list")) GR <- try(as(GR, "GRangesList"))
           if (Sys.info()['sysname'] == "Linux") {
             bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
           } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")

           GR <- bplapply(GR, getGRegionsStat, win.size, step.size, grfeatures,
                   stat, absolute, select.strand, column, prob, entropy,
                   maxgap, minoverlap, scaling, logbase,
                   type, ignore.strand, na.rm, BPPARAM=bpparam)
           return(GR)
           }
)

#' @aliases getGRegionsStat, GRangesList-method
#' @rdname getGRegionsStat-methods
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
setMethod("getGRegionsStat", signature(GR="GRangesList"),
          function(GR, win.size=350, step.size=350, grfeatures=NULL,
                   stat=c("sum", "mean", "gmaean", "median", "density"),
                   absolute=FALSE, select.strand=NULL, column=1L,
                   prob=FALSE, entropy=FALSE, maxgap=-1L, minoverlap=0L,
                   scaling=1000L, logbase = 2,
                   type=c("any", "start", "end", "within", "equal"),
                   ignore.strand=FALSE, na.rm=TRUE, num.cores = 1L, tasks = 0){
            if (Sys.info()['sysname'] == "Linux") {
              bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
            } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")

            GR <- bplapply(GR, getGRegionsStat, win.size, step.size, grfeatures,
                         stat, absolute, select.strand, column, prob, entropy,
                         maxgap, minoverlap, scaling, logbase, type,
                         ignore.strand, na.rm, BPPARAM=bpparam)
            return(GR)
          }
)

