#' @rdname findCutpoint
#' @title Find a cutoff of divergences of methylation level values
#' @description A function to help on the decision of which is the best cutoff
#'     value for DIMP/DMP predictions. The genome-wide methylation changes that
#'     occurs in any living organism is the result of the superposition of
#'     several stochastic processes: the inherent stochasticity of biological
#'     processes and, particular, ultimately, it derives from the stochasticity
#'     of biochemical processes. On this scenario, there is not way to say with
#'     absolute determinism where a given value of an information divergence is
#'     a true positive value or a true negative value. All what we can do is
#'     the estimation of performance indicators like accuracy, sensitivity,
#'     false positive rate, etc., to evaluate the consequences of our decision
#'     on what we consider a true positive or a true negative. For example, a
#'     difference of methylation levels of 100% observed in a treatment group of
#'     samples in given cytosine position does not means that this difference
#'     will not be observed in some sample from the control group. Without any
#'     doubt about it, such a different can be found in control samples as well.
#'     The fluctuation theorem guaranty such an outcome, which in the current
#'     context is a consequence of the action of second law of thermodynamics on
#'     living organisms.
#' @details Given a numerical vector of cutoff values for the divergences of
#'     methylation level values, or p-values cutoffs, this function search for
#'     the cutoff value that yield the best classification performance for the
#'     specified performance indicator.
#' @param LR A list of GRanges, a GRangesList, a CompressedGRangesList object.
#'     Each GRanges object from the list must have at least two columns: a
#'     column containing the total variation of methylation level (TV,
#'     difference of methylation levels) and a column containing a divergence of
#'     methylation levels (it could be TV or  Hellinger divergence) or a
#'     column with a p-value from where the cutpoint will be found
#'     (see example).
#'
#' @param min.tv Minimum value for the total variation distance (TVD; absolute
#'     value of methylation levels differences, \eqn{TVD = abs(TV)}). Only
#'     sites/ranges k with \eqn{TVD_{k} > min.tv} are analyzed. Defaul
#'     min.tv = 0.25.
#' @param tv.cut A cutoff for the total variation distance to be applied to each
#'     site/range. If tv.cut is provided, then sites/ranges k with
#'     \eqn{TVD_{k} < tv.cut} are considered TRUE negatives and
#'     \eqn{TVD_{k} > tv.cut} TRUE postives. Its value must be NULLor a number
#'     \eqn{0 < tv.cut < 1}.
#' @param predcuts A numerical vector of possible cutoff values (cutpoints) for
#'     a divergence of methylation levels value or a p-value, according with the
#'     magnitude given in div.col or in pval.col, respectively. For each
#'     cutpoint k the values greater than predcuts[k] are predicted TRUE
#'     (positives), otherwise are predicted FALSE (negatives).
#' @param tv.col Column number where the total variation is located in the
#'     metadata from each GRanges object.
#' @param div.col Column number for divergence of methylation levels used in the
#'     estimation of the cutpoints. Default: NULL. One of the parameter values
#'     div.col or pval.col must be given.
#' @param pval.col Column number for p-value used in the estimation of the
#'     cutpoints. Default: NULL. One of the parameter values div.col or pval.col
#'     must be given.
#' @param stat An integer number indicating the statistic to be used in the
#'     testing. The mapping for statistic names are:
#'     0 = "All" 1 = "Accuracy", 2 = "Sensitivity", 3 = "Specificity",
#'     4 = "Pos Pred Value", 5 = "Neg Pred Value", 6 = "Precision",
#'     7 = "Recall", 8 = "F1",  9 = "Prevalence", 10 = "Detection Rate",
#'     11 = "Detection Prevalence", 12 = "Balanced Accuracy".
#' @param maximize Whether to maximize the performance indicator given in
#'     parameter 'stat'. Default: TRUE.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see 'bplapply' function from
#'     BiocParallel package).
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom S4Vectors mcols
#' @importFrom caret confusionMatrix
#' @return A list with the classification repformance results for the best
#'     cutoff value in the ranges of predcuts supplied.
#' @export
#' @author Robersy Sanchez
#' @examples
#' # load simulated data of potential methylated signal
#' data(sim_ps)
#'
#' # Vector of cutoff values
#' cuts = c(2, 5, 10, 15, 18, 20, 21, 22, 25, 27, 30, 35, 40,
#'         45, 50, 55, 60)
#
#' # === To find the cutpoint that maximize the accuracy ===
#' pre.cut.acc = findCutpoint(LR = PS, min.tv = 0.25, tv.cut = 0.5,
#'                             predcuts = cuts, tv.col = 7L, div.col = 9,
#'                             stat = 1, num.cores = 15)

findCutpoint <- function(LR, min.tv, tv.cut, predcuts, tv.col, div.col=NULL,
                          pval.col=NULL, stat = 1, maximize = TRUE,
                          num.cores = 1L) {
  # statistic names
  # 0 = "All" 1 = "Accuracy", 2 = "Sensitivity", 3 = "Specificity",
  # 4 = "Pos Pred Value", 5 = "Neg Pred Value", 6 = "Precision", 7 = "Recall",
  # 8 = "F1",  9 = "Prevalence", 10 = "Detection Rate",
  # 11 = "Detection Prevalence", 12 = "Balanced Accuracy"

  require(BiocParallel)
  bpparam <- MulticoreParam(workers=num.cores, tasks=0)
  res <- bplapply(predcuts, function(k)
    classPerform(LR=LR, min.tv=min.tv, tv.cut=tv.cut, cutoff=k, tv.col=tv.col,
                 div.col=div.col, pval.col=pval.col, stat=stat),
       BPPARAM = bpparam)
  res <- unlist(res)
  if (!all(is.na(res))) {
    if (maximize) cutp = predcuts[which.max(res)] else
      cutp = predcuts[which.min(res)]
    perf <- classPerform(LR=LR, min.tv=min.tv, tv.cut=tv.cut, cutoff=cutp,
                         tv.col=tv.col, div.col=div.col, pval.col=pval.col,
                         stat=0)

    res = list(optimCutpoint=cutp, Statistic=res, Performance=perf[[1]],
               FDR=perf[[2]])
  } else res = list(optimCutpoint=NA, Statistic=NA, Performance=NA,
                    FDR=NA)
  return(res)
}

