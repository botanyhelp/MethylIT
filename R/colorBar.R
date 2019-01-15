#' @rdname colorBar
#' @title Color bar for heatmaps
#' @description This function is for internal use. It creates a color scaled
#'     bar for heatmap.
#' @param z matrix of values used to visualize the heatmap
#' @param zlim Limit for the numerical(color) scale, which must be consistent
#'     with parameter 'break'.
#' @param col Palette of color to use
#' @param breaks Numerical vector with the breaks used to produce the heatmap
#' @param horiz Whether to the color bar will be horizontal(= TRUE) or
#'     vertical(= FALSE)
#' @param ylim User defined limits for y-axis. Depending on the orientation, x-
#'     or y-limits may be defined that are different from the z-limits and will
#'     reduce the range of colors displayed.
#' @param xlim User defined limits for x-axis
#' @param lwd Line width.
#' @param cex.axis Cex values for color bar axis labels.
#' @param ... Additional parameter to pass to 'par' R function
#' @return Image with color bar
#'
.colorBar <- function(z, zlim, col=heat.colors(12), breaks, horiz=TRUE,
                       ylim=NULL, xlim=NULL, lwd=0.5, cex.axis=1, ...){

   if (!missing(breaks)) {
       if (length(breaks) != (length(col) + 1)) {
           stop("must have one more break than colour")
       }
   }
   if (missing(breaks) & !missing(zlim)) {
       breaks <- seq(zlim[1], zlim[2], length.out=(length(col) + 1))
   }
   if (missing(breaks) & missing(zlim)) {
       zlim <- range(z, na.rm=TRUE)
       ## adds a bit to the range in both directions
       zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (1E-3)
       zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (1E-3)
       breaks <- seq(zlim[1], zlim[2], length.out=(length(col) + 1))
   }

   poly <- vector(mode="list", length(col))
   for (i in seq(poly)) {
       poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1], breaks[i])
   }

   if (horiz) {YLIM <- c(0,1); XLIM <- range(breaks)}
   if (!horiz) {YLIM <- range(breaks); XLIM <- c(0, 1)}
   if (missing(xlim)) xlim <- XLIM
   if (missing(ylim)) ylim <- YLIM

   plot(x=0.5, y=0.5, type="n", ylim=ylim, xlim=xlim, xaxt="n",
           yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",
           lwd=lwd, ...)
   if (!horiz) axis(side=2, ylim=ylim, xlim=xlim, lwd=lwd, cex.axis=cex.axis, ... )
   if (horiz) axis(side=1, ylim=ylim, xlim=xlim, lwd=lwd, cex.axis=cex.axis, ... )

   for (i in seq(poly)) {
       if (horiz) {
           polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
       }
       if (!horiz) {
           polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
       }
   }
}

