#' @rdname classPerform
#'
#' @title Classification performance based on divergences of methylation levels
#' @description The classification performance based on an information
#'     divergence  (e.g., Hellinger divergence) carried in a list of GRanges
#'     objects. The total variation distance (TVD, absolute difference of
#'     methylation levels) is used as pivot to specify the cytosine sites
#'     considered as true positives and true negatives. Function
#'     \code{\link[caret]{confusionMatrix}} from package "caret" is applied to
#'     get the classification performance.
#' @details Samples from each group are pooled according to the statistic
#'     selected (see parameter pooling.stat) and a unique GRanges object is
#'     created with the methylated and unmathylated read counts for each group
#'     (control and treatment) in the metacolumn. So, a contingence table can be
#'     built for range from GRanges object.
#' @param LR A list of GRanges, a GRangesList, a CompressedGRangesList object.
#'     Each GRanges object from the list must have two columns: methylated
#'     (mC) and unmethylated (uC) counts. The name of each element from the
#'     list must coincide with a control or a treatment name.
#' @param min.tv Minimum value for the total variation distance (TVD; absolute
#'     value of methylation levels differences, \eqn{TVD = abs(TV)}).
#'     Only sites/ranges k with \eqn{TVD_{k} > min.tv} are analyzed. Defaul
#'     min.tv = 0.25.
#' @param tv.cut A cutoff for the total variation distance to be applied to each
#'     site/range. If tv.cut is provided, then sites/ranges k with
#'     \eqn{TVD_{k} < tv.cut} are considered TRUE negatives and
#'     \eqn{TVD_{k} > tv.cut} TRUE postives. Its value must be NULLor a number
#'     \eqn{0 < tv.cut < 1}.
#' @param cutoff A divergence of methylation levels or a p-value cutoff-value
#'     for the the magnitude given in div.col or in pval.col, respectively
#'     (see below). The values greater than 'cutoff' are predicted TRUE
#'     (positives), otherwise are predicted FALSE (negatives).
#' @param div.col Column number for divergence variable used in the performance
#'     analysis and estimation of the cutpoints. Default: NULL. One of the
#'     parameter values div.col or pval.col must be given.
#' @param pval.col Column number for p-value used in the performance
#'     analysis and estimation of the cutpoints. Default: NULL. One of the
#'     parameter values div.col or pval.col must be given.
#' @param stat An integer number indicating the statistic to be used in the
#'     testing. The mapping for statistic names are:
#'     0 = "All" 1 = "Accuracy", 2 = "Sensitivity", 3 = "Specificity",
#'     4 = "Pos Pred Value", 5 = "Neg Pred Value", 6 = "Precision",
#'     7 = "Recall", 8 = "F1",  9 = "Prevalence", 10 = "Detection Rate",
#'     11 = "Detection Prevalence", 12 = "Balanced Accuracy".
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom S4Vectors mcols
#' @importFrom caret confusionMatrix
#' @return A list with the classification repformance results
#' @export
#' @author Robersy Sanchez
#' @examples
#' # load simulated data of potential methylated signal
#' data(sim_ps)
#'
#' classPerform(LR = PS, min.tv = 0.25, tv.cut = 0.4,
#'              cutoff = 68.7, tv.col = 7L, div.col = 9, stat = 0)
classPerform <- function(LR, min.tv = 0.25, tv.cut, cutoff, tv.col,
                         div.col=NULL, pval.col=NULL, stat=1) {
  if (is.null(div.col) && is.null(pval.col))
    stop("*** One of the parameters 'div.col' or 'pval.col' must be not NULL")
  if (is.null(div.col)) target.col = pval.col else target.col = div.col
  if (inherits(LR, "list")) LR <- try(as(LR, "GRangesList"), silent=TRUE)
  if (inherits(LR, "try-error")) stop("*** LR is not a list of GRanges objects")
  if (inherits(LR, "GRangesList")) LR <- unlist(LR)

  LR <- LR[abs(mcols(LR[, tv.col])[,1]) > min.tv]
  ref.signal <- rep("No Event", length(LR))
  idx = which(abs(mcols(LR[, tv.col])[,1]) > tv.cut)
  ref.signal[idx] <- "Event"

  pred.signal <- rep("No Event", length(LR))
  if (is.null(div.col)) {
    idx = which(abs(mcols(LR[, target.col])[,1]) < cutoff)
  } else idx = which(abs(mcols(LR[, target.col])[,1]) > cutoff)
  pred.signal[idx] <- "Event"

  tbl <- table(pred.signal, ref.signal)

  if (sum(dim(tbl)) > 2) {
    if (!is.element("No Event", unique(pred.signal))) {
      tbl <- rbind(tbl, c(0,0))
      rownames(tbl) <- c("Event", "No Event")
      attributes(tbl) <- list(dim=c(2,2),
                              dimnames=list(pred.signal=c("Event", "No Event"),
                                            ref.signal=c("Event", "No Event")),
                              class="table")
    }

    if (!is.element("Event", unique(pred.signal))) {
      tbl <- rbind(c(0,0), tbl)
      rownames(tbl) <- c("Event", "No Event")
      attributes(tbl) <- list(dim=c(2,2),
                              dimnames=list(pred.signal=c("Event", "No Event"),
                                            ref.signal=c("Event", "No Event")),
                              class="table")
    }
  }

  if (sum(dim(tbl) == 2)) {
    conf.mat <- try(confusionMatrix(tbl, positive = "Event"), silent=TRUE)
    if (inherits(conf.mat, "try-error"))
       stop("*** It was no possible to build a 2x2 confusion matrix", "\n",
            "Perhaps parameters 'tv.cut' and 'cutoff' have too high values")
    if (stat > 1) {
      res <- conf.mat$byClass[stat - 1]
    } else {
      if (stat == 1) res <- conf.mat$overall[1]
      if (stat == 0) {
        FDR = tbl[1,2]/sum(tbl[1,])
        res <- list(Performance = conf.mat, FDR = FDR)
      }
    }
  } else res <- NA

  return(res)
}
