#' @rdname dmpDensity
#' @title Linear density of DMPs at a given genomic region
#' @description The linear density of DMPs in a given genomic region (GR) is
#'     defined according with the classical terminology in physics, i.e., as the
#'     measure of the physical quantity of any characteristic value per unit of
#'     length. In the current case, as the amount of DIMPs per nucleotide base.
#'
#' @details Since the number of DIMPs along the DNA sequence vary, the local
#'     density of DMPs \eqn{\rho_i} at a fixed interval \eqn{\Delta} l_i is
#'     defined by the quotient \eqn{\rho_i = \Delta DMP_i/\Delta l_i} is the
#'     amount of DIMPs at the fixed interval. Likewise the local density of
#'     non-DIMPs is defined as \eqn{\rho_i = \Delta nonDMP_i/\Delta l_i}.
#'     Notice that for a specified methylation context, e.g., CG,
#'     \eqn{\Delta CG_i - \Delta DMP_i}, where \eqn{\Delta CG} is the
#'     amount CG positions at the given interval. The linear densities are
#'     normalized as \eqn{\rho_i/\rho_max}, where \eqn{\rho_max} is the maximum
#'     of linear density found in a given GR.
#' @param GR A genomic GRanges object carrying the genomic region where the
#'     estimation of the DMP density will be accomplished.
#' @param cut.col Integer denoting the GR metacolumn where the decision variable
#'     about whether a position is DMP is located. Default cut.col = 1.
#' @param cutoff Cut value to decide wheter the value of the variable used to
#'     estimate the density is a DMP at each position. If missing, then
#'     cutoff is estimated as the first queantile greater than zero from the
#'     values given in the GR column \emph{cut.col}.
#' @param Chr A character string. Default NULL. If the GR object comprises
#'     several chromosomes, then one chromosome must be specified. Otherwise the
#'     density of first chromosome will be returned.
#' @param start.pos,end.pos Start and end positions, respectively, of the GR
#'     where the density of DMPs will be estimated. Default NULL. If NULL
#'     densities will be estimated for the whole GR and the specified
#'     chromosome.
#' @param int.size1,int.size2 The interval/window size where the density of DMP
#'     and no DMPs are computed. Default Null.
#' @param breaks Integer. Number of windows/intervals to split the GR. Deafult
#'     NULL. If provided, then it is applied to compute the densities of DMPs
#'     and no-DMPs. If 'int.size1', 'int.size2', and 'breaks' are NULL, then the
#'     breaks are computed as:
#'     \code{breaks <- min(150, max(start(x))/nclass.FD(start(x)),
#'     na.rm = TRUE)},
#'     where function \emph{nclass.FD} (\code{\link[grDevices]{nclass}}) applies
#'     Freedman-Diaconis algorithm.
#' @param scaling Logic value to deside whether to perform the scaling of the
#'     estimated density values or not. Default is TRUE.
#' @param plot Logic. Whether to produce a grahic or not. Default, plot = TRUE.
#' @param noDMP.dens Logic whether to produce the graphics for no-DMP density.
#'     Default is TRUE
#' @param xlabel X-axis label. Default \emph{xlabel = "Coordinate"}.
#' @param ylabel Y-axis label. Default \emph{ylabel = "Normalized density"}.
#' @param col.dmp Color for the density of DMPs in the graphic.
#' @param col.ndmp Color for the density of no DMPs in the graphic.
#' @param yintercept If plot == TRUE, this is the position for an horizantal
#'     line that intercept the y-axis. Default yintercept = 0.25.
#' @param col.yintercept Color for the horizantal line 'yintercept'. Default
#'     \emph{col.yintercept = 'blue'}
#' @param type.yintercept Line type for the horizantal line 'yintercept'.
#'     Default \emph{type.yintercept = "dashed"}.
#' @param dig.lab integer which is used when labels are not given. It determines
#'     the number of digits used in formatting the break numbers.
#' @examples
#' set.seed(349)
#' ## An auxiliary function to generate simulated hypothetical values from a
#' ## variable with normal distribution
#' hypDT <- function(mean, sd, n, num.pos, noise) {
#'     h <- hist(rnorm(n, mean = mean, sd = sd), breaks = num.pos, plot = FALSE)
#'     hyp <- h$density * 60 + runif(length(h$density)) * noise
#'     return(hyp)
#' }
#'
#' ## To generate a matrix of values with variations introduced by noise
#' hyp <- hypDT(mean = 5, sd = 30, n = 10^5, noise = 4, num.pos = 8000)
#' ## A GRanges object is built, which will carries the previous matrix on its
#' ## meta-columns
#' l <- length(hyp)
#' starts <- seq(0, 30000, 3)[1:l]
#' ends <- starts
#' GR <- GRanges(seqnames = "chr1", ranges = IRanges(start = starts,
#'                 end = ends))
#' mcols(GR) <- data.frame(signal = hyp)
#'
#' # If plot is TRUE a grphic is printed. Otherwise data frame is returned.
#' p <- dmpDensity(GR, plot = FALSE)
#'
#' # If ggplot2 package is installed, then graphic can customized using
#' # the returned data frame 'p':
#'
#' # library(ggplot2)
#' ## Auxiliar function to write scientific notation in the graphics
#' # fancy_scientific <- function(l) {
#' #   #'turn in to character string in scientific notation
#' #   l <- format( l, scientific = TRUE, digits = 1 )
#' #   l <- gsub("0e\\+00","0",l)
#' #   #'quote the part before the exponent to keep all the digits
#' #   l <- gsub("^(.*)e", "'\\1'e", l)
#' #   #'turn the 'e+' into plotmath format
#' #   l <- gsub("e", "%*%10^", l)
#' #   l <- gsub("[+]", "", l )
#' #   #'return this as an expression
#' #   parse(text=l)
#' # }
#' #
#' # max.pos = max(p$DMP.coordinate)
#' # ggplot(data=p) +
#' #   geom_line(aes(x=DMP.coordinate, y=DMP.density), color="red") +
#' #   geom_hline(aes(yintercept=0.25), linetype="dashed",
#' #              colour="blue", show.legend=FALSE ) +
#' #   geom_line(aes(x=coordinate, y=density), color="blue") +
#' #   xlab("Coordinate") + ylab("Normalized density") +
#' #   scale_y_continuous(breaks=c(0.00, 0.25, 0.50, 0.75, 1.00)) +
#' #   scale_x_continuous(breaks=c(0.00, 0.25 *max.pos, 0.50*max.pos,
#' #                               0.75*max.pos, max.pos),
#' #                      labels = fancy_scientific) +
#' #   expand_limits(y=0)
#' @importFrom GenomicRanges setdiff
#' @importFrom grDevices nclass.FD
#' @return If plot is TRUE will return a graphic with the densities of DMPs and
#'     and no DMPs. If plot is FALSE a data frame object with the density of
#'     DMPs and not DMPs will be returned.
#' @author Robersy Sanchez
#' @export

dmpDensity <- function(GR, column=1, cut.col=1, cutoff, Chr=NULL,
                       start.pos=NULL, end.pos=NULL, int.size1=NULL,
                       int.size2=NULL, breaks=NULL, scaling=TRUE, plot=FALSE,
                       noDMP.dens=TRUE, xlabel="Coordinate",
                       ylabel="Normalized density", col.dmp="red",
                       col.ndmp="blue", yintercept=0.25,
                       col.yintercept="magenta", type.yintercept="dashed",
                       dig.lab=3) {
   if (class(GR) != "GRanges") stop("GR must be a 'GRanges' object")

   # Auxiliar function to write scientific notation in the graphics
   fancy_scientific <- function(l) {
       # turn in to character string in scientific notation
       l <- format( l, scientific = TRUE, digits = 1 )
       l <- gsub("0e\\+00","0",l)
       # quote the part before the exponent to keep all the digits
       l <- gsub("^(.*)e", "'\\1'e", l)
       # turn the 'e+' into plotmath format
       l <- gsub("e", "%*%10^", l)
       l <- gsub("[+]", "", l )
       # return this as an expression
      parse(text=l)
   }
   # Function to compute the density of DMPs
   points <- function(x, int.length=NULL, scale, dig.lab, breaks=NULL) {
       intv = start(x)
       # Split the gene into fixex interval
       if (is.null(int.length) & is.null(breaks)) {
           breaks <- min(150, max(start(x))/nclass.FD(start(x)), na.rm = TRUE)
       }
       if (!is.null(int.length)) breaks <- round(length(x)/int.length)
       dmp <- table(cut(intv, breaks=breaks, dig.lab=dig.lab))
       # Middle interval position
       pos <- rowMeans(t(unname(
               sapply(gsub("[(]*[]]*", "", names(dmp)),
                      function(s) as.numeric(strsplit(s, split = "," )[[1]])))))
       dmp <- as.numeric(dmp)
       if(scaling) dmp <- dmp/max(dmp, na.rm = TRUE)
       return(data.frame(pos, dmp))
   }

   if (!is.null(Chr)) GR <- GR[ seqnames(GR) == Chr,] else {
       chr <- unique(seqnames(GR))
       if(length(chr) > 1) GR <- GR[ seqnames(GR) == chr[1],]
   }

   if (!is.null(start.pos) && !is.null(end.pos)) {
       starts = start(GR)
       if (end.pos > max(starts))
           stop("'end.pos' is greater than the higest coordinate in the GR")
       if (end.pos > max(starts))
           stop("'start.pos' is greater than the lowest coordinate in the GR")
       idx1 = min(which(starts >= start.pos))
       idx2 = max(which(start(GR) <= end.pos))
       GR <- GR[idx1:idx2]
   }

   if (missing(cutoff)) {
     qs <- summary(mcols(GR)[,1])[2:6]
     cutoff <- qs[qs > 0][1]
   }

   idx <- (mcols(GR)[, cut.col] >= cutoff)
   x1 <- GR[idx, ]
   x2 <- GenomicRanges::setdiff(GR, x1, ignore.strand=TRUE)
   x1 <- points(x=x1, int.length=int.size1, scale=scale, dig.lab=dig.lab,
                breaks=breaks)
   x2 <- points(x=x2, int.length=int.size2, scale=scale, dig.lab=dig.lab,
                breaks=breaks)
   if (plot) {
         max.pos = max(x1$pos)
         tick.pos=c(0.00, 0.25 *max.pos, 0.50*max.pos, 0.75*max.pos, max.pos)
         par(las=1)
         plot(x1$pos, x1$dmp, type="l", xlab = xlabel, xaxt="n",
              ylab=ylabel, col=col.dmp, ylim=c(0,1))
         axis(1, at=tick.pos, label=fancy_scientific(tick.pos))
         if (noDMP.dens) lines(x2$pos, x2$dmp, col=col.ndmp)
         abline(h=yintercept, col=col.yintercept, lty = type.yintercept)
   } else {
     return(data.frame(DMP.coordinate=x1$pos, DMP.density=x1$dmp,
                       coordinate=x2$pos, density=x2$dmp))
   }
}


