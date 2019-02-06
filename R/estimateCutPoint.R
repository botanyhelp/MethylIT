#' @rdname estimateCutPoint
#'
#' @title Estimate cutpoints to distinguish the treatment methylation signal
#'     from the control
#' @description Given a list of two GRanges objects, control and treatment,
#'     carrying the potential signals (prior classification) from controls and
#'     treatments in terms of an information divergence
#'     (given the meta-columns), the function estimates the cutpoints of the
#'     control group versus treatment group.
#' @details The function performs an estimation of the optimal cutpoint from
#'     the area under the curve (AUC) of a receiver operating characteristic
#'     (ROC) and the accuracy for the classification of the cytosine positions
#'     based on the cutpoint into two classes: DIMPs from control and DIMPs
#'     from treatment.
#'
#'     In this context, the AUC is the probability of being able to distinguish
#'     the biological regulatory signal naturally generated in the control from
#'     that one induced by the treatment. The cytosine sites carrying a
#'     methylation signal shall be called differentially informative methylated
#'     positions (DIMP). Now, the probability that a DIMP is not induced by the
#'     treatment is given by the probability of false alarm
#'     (PFA, false positive). That is, the biological signal is naturally
#'     present in the control as well as in the treatment.
#'
#' @param LR A list of GRanges objects containing a divergence variable used to
#'     perform ROC analysis and estimate the cutpoint
#' @param control.names Names/IDs of the control samples. Each GRanges object
#'     must correspond to a sample, for example, sample 's1'. Then this sample
#'     can be accessed in the list of GRanges objects as LR$s1.
#' @param treatment.names Same type and function as 'control.names'.
#' @param div.col Column number for divergence variable used in the ROC
#'     analysis and estimation of the cutpoints.
#' @param absolute Logic (default, FALSE). Total variation (TV, the difference
#'     of methylation levels) is normally an output in the downstream MathylIT
#'     analysis. If 'absolute = TRUE', then TV is tranformed into |TV|, which is
#'     an information divergence that can be fitted to Weibull or to Generalized
#'     Gamma distribution. So, if the nonlinear fit was performed for |TV|, then
#'     here absolute must be turned TRUE.
#' @param grouping Logic (default, FALSE). If TRUE, then all control samples are
#'     grouped into one set, and the same for the treatment samples. The new
#'     samples are named: 'ctrl' and 'treat' and, consequently, only a cutpoint
#'     (treat' versus 'ctrl') is estimated.
#'@param verbose If TRUE, prints the function log to stdout.
#' @return A list of three matrices cutpoint matrix values, AUC matrix values,
#'     and accuracy matrices values. These matrices values derives from all
#'     possible ROC analysis: control sample_i versus treatment sample_j (i,j =
#'     1,2, ...).
#'
#' @examples
#'
#'     set.seed(123)
#'     num.points <- 1000
#'
#'     ## A list of GRanges objects with simulated Hellinger divergences in
#'     ## their metacolumns.
#'     HD <- GRangesList(
#'         sample1 = makeGRangesFromDataFrame(
#'                     data.frame(chr = "chr1", start = 1:num.points,
#'                             end = 1:num.points, strand = '*',
#'                             hdiv = rweibull(1:num.points, shape = 0.45,
#'                                     scale = 0.2)),
#'                     keep.extra.columns = TRUE),
#'         sample2 = makeGRangesFromDataFrame(
#'                     data.frame(chr = "chr1", start = 1:num.points,
#'                             end = 1:num.points, strand = '*',
#'                             hdiv = rweibull(1:num.points, shape = 0.75,
#'                                     scale = 1)),
#'                     keep.extra.columns = TRUE))
#'
#'     ## Nonlinear fit of Weiblul distribution
#'     nlms <- nonlinearFitDist(HD, column = 1, verbose = FALSE)
#'
#'     ## Estimation of the potential signal and cutpoints
#'     PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 1, alpha = 0.05)
#'     cutpoints <- estimateCutPoint(PS, control.names = "sample1",
#'                                     treatment.names = c("sample2"),
#'                                     div.col = 1, verbose = FALSE)
#' # Let's add an empty GRanges
#' PS$sample1.1 <- GRanges()
#' # The empty GRanges will be ignored
#' cutpoints <- estimateCutPoint(PS, control.names = c("sample1", "sample1.1"),
#'                               treatment.names = c("sample2"),
#'                               div.col = 1, verbose = FALSE)
#' # Let's make both control GRanges empty.
#' PS$sample1 <- GRanges()
#' # Warning is given
#' cutpoints <- estimateCutPoint(PS, control.names = c("sample1", "sample1.1"),
#'                               treatment.names = c("sample2"),
#'                               div.col = 1, verbose = FALSE)
#' # Let's make the treatment GRange empty.
#' PS$sample2 <- GRanges()
#' # The error is reported
#' cutpoints <- estimateCutPoint(PS, control.names = c("sample1", "sample1.1"),
#'                               treatment.names = c("sample2"),
#'                               div.col = 1, verbose = FALSE)
#' @export
estimateCutPoint <- function(LR, control.names, treatment.names, div.col=NULL,
                             absolute=FALSE, grouping=FALSE, verbose = TRUE) {

   sn <- names(LR)

   infDiv <- function(LR) {
       ## This builds data frames from the list or ranges LR
       ## to be used for ROC analysis
       ## LR: list of sample GRanges
       if (is.null(div.col))
           stop(paste("* Provide the column number for divergence variable ",
                   "('div.col')"))
       sn <- names(LR)
       dt <- list()
       for (k in 1:length(LR)) {
           x = LR[[k]]@elementMetadata[ , div.col]
           if (absolute) x = abs(x)
           dt[[k]] <- data.frame(idiv=x, treat=sn[k])
       }
       names( dt ) <- sn
       return(dt)
   }

   roc <- function(dt0, dt1) {
       ## This function build the ROC and estimate the
       ## Hellinger divergence cutoff point starting from
       ## which TRUE DIMPs are found.
       ## dt0 & dt1: data frames built with HD function
       ## dt0: control
       ## dt1: treatment
       ## folder: place where the ROC will be printed
       dt <- rbind(dt0, dt1)
       l <- levels(dt$treat)
       dt$status <- as.character(dt$treat)
       dt$status[dt$status == l[1]] <- 0
       dt$status[dt$status == l[2]] <- 1
       dt$status <- as.numeric(dt$status)
       rownames(dt) <-NULL

       m <- as.matrix(base::table(dt$idiv,  dt$status))
       m <- addmargins(rbind(0, m), 2)
       nr <- nrow(m)
       m <- apply(m, 2, cumsum)
       sens <- (m[nr, 2] - m[, 2])/m[nr, 2]
       spec <- m[, 1]/m[nr, 1]
       # res <- data.frame(sens, spec, prdtr = as.numeric(rownames(m)))
       auc <- sum((sens[-1] + sens[-nr])/2 * abs(diff(1 - spec)))
       idx <- which.max(sens[-1] + spec[-nr])
       cutpoint <- as.numeric(names(idx))
       conf.matrix <- table(dt$idiv > cutpoint, dt$status)
       list(cutpoint=cutpoint, acc=sum(diag(conf.matrix)) / sum(conf.matrix),
           AUC=auc)
   }

   lcc <- unlist(lapply(control.names, function(k) length(LR[[k]]) > 0))
   ltt <- unlist(lapply(treatment.names, function(k) length(LR[[k]]) > 0))


   if (sum(ltt) < length(LR[treatment.names])) {
       if (sum(ltt) == 0) {
           text <- paste0("All the GRanges objects from treatment group are ",
                       "empty, please check your data")
           stop(text)
       } else treatment.names <- treatment.names[ltt]
   }

   if (sum(lcc) == 0) {
       LR <- LR[treatment.names]
       min.div <- min(unlist(lapply(LR, function(l)
           min(l@elementMetadata[, div.col], na.rm = TRUE))), na.rm = TRUE)
       min.div <- min.div[min.div > 0]
       text <- c("All the GRanges objects from your control are empty ", "\n",
                       "So, the cutpoint is the min(div) > 0 value found", "\n",
                       "in the treatment group")
       warning(text)
       return(list(cutpoint=min.div, auc=NA, accuracy=NA))
   } else {
       control.names <- control.names[lcc]
       LR <- LR[c(control.names, treatment.names)]
       if (grouping) {
           LR = list(unlist(as(LR[control.names], "GRangesList")),
                   unlist(as(LR[treatment.names], "GRangesList")))
           sn <-c("ctrl", "treat")
           names(LR) <- sn
           control.names <- "ctrl"
           treatment.names <- "treat"
       }
       LR <- infDiv(LR)
       idx.ct <- match(control.names, sn)
       idx.tt <- match(treatment.names, sn)

       cutpoint <- c()
       accuracy <- c()
       auc <- c()
       for (j in idx.ct) {
           for (k in idx.tt) {
               if (verbose)
                   message("*** Comparison: ", paste0(sn[j], " versus ", sn[k]))
               res <- roc(dt0=LR[[j]], dt1=LR[[k]])
               cutpoint <- c(cutpoint, res$cutpoint)
               accuracy <- c(accuracy, res$acc)
               auc <- c(auc, res$AUC)
           }
       }

       accuracy <- data.frame(matrix(accuracy, nrow=length(idx.tt)))
       cutpoint <- data.frame(matrix(cutpoint, nrow=length(idx.tt)))
       auc <- data.frame(matrix(auc, nrow=length(idx.tt)))

       colnames(accuracy) <- control.names
       rownames(accuracy) <- treatment.names
       colnames(cutpoint) <- control.names
       rownames(cutpoint) <- treatment.names
       colnames(auc) <- control.names
       rownames(auc) <- treatment.names

       return(list(cutpoint=cutpoint, auc=auc, accuracy=accuracy))
   }
}
