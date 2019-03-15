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
#' @importFrom S4Vectors mcols
#' @export
estimateCutPoint <- function(LR, control.names, treatment.names, simple = TRUE,
                       column=c(hdiv=TRUE, TV=TRUE, wprob=FALSE, pos=FALSE),
                       classifier1=c("logistic", "pca.logistic", "lda",
                                   "qda","pca.lda", "pca.qda"),
                       classifier2=NULL, tv.cut = 0.25, div.col = NULL,
                       post.cut = 0.5, clas.perf = FALSE, prop=0.6, n.pc=1,
                       find.cut=FALSE, cut.interval = c(0.5, 0.8),
                       cut.incr = 0.01, stat = 1, num.cores=1L, tasks=0L,
                       tol = .Machine$double.eps^0.5, verbose = TRUE, ...) {
   if (!simple && sum(column) == 0) {
       cat("\n")
       stop(paste("*** At least one of columns with the predictor \n",
               "variables: 'hdiv', 'TV', 'wprob', or 'pos' must be provided"))
   }
   if ((classifier1[1] != "logistic" ) && sum(column) < n.pc) {
       cat("\n")
       stop(paste("* The number of predictor variables must be greater \n",
               "or equal to n.pc"))
   }

   if (is.null(classifier2[1])) {
       if ((classifier2[1] != "logistic" ) && sum(column) < n.pc) {
           stop(paste("* The number of predictor variables must be greater or ",
                       "equal to n.pc"))
       }
   }

   # --------------------------- valid "pDMP" object-------------------------- #
   validateClass(LR)
   # ------------------------------------------------------------------------- #
   sn <- names(LR)

   # -------------------------- Auxiliary functions -------------------------- #
   infDiv <- function(LR, div.col = NULL) {
       ## This builds data frames from the list or ranges LR
       ## to be used for ROC analysis
       ## LR: list of sample GRanges
       if (is.null(div.col))
           stop(paste("* Provide a divergence column"))

       sn <- c("ctrl", "treat")
       dt <- list()
       for (k in 1:2) {
           x = LR[[k]]@elementMetadata[, div.col]
           dt[[k]] <- data.frame(idiv = abs(x), treat = sn[k])
       }
       names(dt) <- sn
       return(dt)
   }
   roc <- function(dt) {
       ## This function build the ROC and estimate the
       ## Hellinger divergence cutoff point starting from
       ## which TRUE DIMPs are found.
       ## dt0 & dt1: data frames built with HD function
       ## dt0: control
       ## dt1: treatment
       ## folder: place where the ROC will be printed
       dt <- do.call(rbind, dt)
       l <- levels(dt$treat)
       dt$status <- as.character(dt$treat)
       dt$status[dt$status == l[1]] <- 0
       dt$status[dt$status == l[2]] <- 1
       dt$status <- as.numeric(dt$status)
       rownames(dt) <- NULL

       m <- as.matrix(base::table(dt$idiv,  dt$status))
       m <- addmargins(rbind(0, m), 2)
       nr <- nrow(m)
       m <- apply(m, 2, cumsum)
       sens <- (m[nr, 2] - m[, 2])/m[nr, 2]
       spec <- m[, 1]/m[nr, 1]
       auc <- sum((sens[-1] + sens[-nr])/2 * abs(diff(1 - spec)))
       # Youden Index
       idx <- which.max(sens[-1] + spec[-nr])
       cutpoint <- as.numeric(names(idx))
       conf.matrix <- table(dt$idiv > cutpoint, dt$status)
       list(cutpoint=cutpoint, acc=sum(diag(conf.matrix)) / sum(conf.matrix),
           AUC=auc)
   }

   res <- list(cutpoint = NA,
               testSetPerformance = NA,
               testSetModel.FDR = NA,
               model = NA,
               modelConfMatrix = NA,
               initModel = NA,
               postProbCut = NA,
               classifier = NA,
               statistic = NA,
               optStatVal = NA
   )
   res <- structure(res, class = c("CutPoint", "list"))

   # ------------------------------------------------------------------------- #
   if (!is.null(control.names)&&!is.null(treatment.names))
     LR = try(LR[c(control.names, treatment.names)], silent=TRUE)
   if (inherits(LR, "try-error"))
     stop("List's names does not match control & treatment names")

   lcc <- unlist(lapply(control.names, function(k) length(LR[[k]]) > 0))
   ltt <- unlist(lapply(treatment.names, function(k) length(LR[[k]]) > 0))
   vn <- c("hdiv", "TV")
   vn <- vn[match(TRUE, column[vn])]
   if (is.null(div.col)) div.col <- vn

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
       return(cutpoint = min.div)
   } else {
       control.names <- control.names[lcc]
       LR <- LR[c(control.names, treatment.names)]

       # ========================== Grouping ==================================
       LR = list(unlist(LR[control.names]), unlist(LR[treatment.names]))
       names(LR) <- c("ctrl", "treat")
       LR <- structure(LR, class = 'pDMP')
       tv.col = match("TV", colnames(mcols(LR[[1]])))
       classes <- c(rep("CT", length(LR$ctrl)),
                    rep("TT", length(LR$treat)))

       if (simple) {
           reslt <- roc(dt = infDiv(LR, div.col = div.col))
           cutpoint <- reslt$cutpoint
           predClasses <- unlist(LR)$hdiv > cutpoint
           predClasses[ predClasses == TRUE ] <- "TT"
           predClasses[ predClasses == FALSE ] <- "CT"

           if (clas.perf && !find.cut) {
               dmps <- selectDIMP(LR, div.col = div.col, cutpoint = cutpoint,
                                   tv.col = tv.col, tv.cut = tv.cut)

               conf.mat <- evaluateDIMPclass(dmps, column = column,
                                           control.names = "ctrl",
                                           treatment.names = "treat",
                                           classifier=classifier1[1],
                                           prop=prop,
                                           output = "conf.mat",
                                           num.cores=num.cores,
                                           tolerance=tol,
                                           tasks=tasks, verbose = FALSE, ...)

               predClasses <- predict(object = conf.mat$model, newdata = LR)

               res$cutpoint <- cutpoint
               res$testSetPerformance <- conf.mat$Performance
               res$testSetModel.FDR <- conf.mat$FDR
               res$model <- conf.mat$model
               res$modelConfMatrix <- table(classes, predClasses)
               res$initModel <- "Youden Index"
               res$classifier <- classifier1[1]
           } else {
               res$cutpoint <- cutpoint
               res$modelConfMatrix <- table(classes, predClasses)
               res$acc <- reslt$acc
               res$auc <- res$AUC
               res$initModel <- "Youden Index"
           }
       } else {
       # ====================== Fit classifier model 1 ========================
       # ------------------------------------------------------------------- #
       # Auxiliar function to find cutpoint/intersection point of the two
       # gamma distributions
           cutFun <- function(divs, post.cut) {
               idx <- which(post > post.cut)
               TT <- divs[idx]
               CT <- divs[-idx]
               CT <- unlist(CT)
               TT <- unlist(TT)
               return(max(c(max(CT$hdiv),min(TT$hdiv))))
           }
       # --------------------------------------------------------------------- #
           if (!find.cut) {
               if (is.null(classifier2)) classifier2 <- classifier1
               conf.mat <- evaluateDIMPclass(LR, column = column,
                                           control.names = "ctrl",
                                           treatment.names = "treat",
                                           classifier = classifier1[1],
                                           prop = prop,
                                           output = "conf.mat",
                                           num.cores=num.cores,
                                           tolerance=tol,
                                           tasks=tasks, verbose = FALSE, ...)

               post <- predict(object = conf.mat$model, newdata = LR,
                               type = "posteriors")
               cutpoint <- cutFun(divs = unlist(LR), post.cut)
               dmps <- selectDIMP(LR, div.col = div.col, cutpoint = cutpoint,
                                   tv.col = tv.col, tv.cut = tv.cut)
               if (length(dmps$ctrl) == 0) {
                   warning("For the estimated cutpoint = ", cutpoint,
                       ", only treatment's DMPs are found. \n",
                       "A classification model with posterior probability = \n",
                       "0.5 will be applied", sep = "")

                   predClasses <- predict(object = conf.mat$model, newdata = LR)
                   res$postProbCut <- 0.5
                   res$cutpoint <- cutFun(divs = unlist(LR), 0.5)
                   res$testSetPerformance <- conf.mat$Performance
                   res$testSetModel.FDR <- conf.mat$FDR
                   res$model <- conf.mat$model
                   res$modelConfMatrix <- table(classes, predClasses)
                   res$initModel <- classifier1[1]
                   res$classifier <- classifier1[1]
               } else {
                   conf.mat <- evaluateDIMPclass(dmps, column = column,
                                               control.names = "ctrl",
                                               treatment.names = "treat",
                                               classifier = classifier2[1],
                                               prop = prop,
                                               output = "conf.mat",
                                               num.cores=num.cores,
                                               tolerance=tol,
                                               tasks=tasks, verbose = FALSE,
                                               ...)

                   predClasses <- predict(object = conf.mat$model, newdata = LR)

                   res$cutpoint <- cutpoint
                   res$postProbCut <- post.cut
                   res$testSetPerformance <- conf.mat$Performance
                   res$testSetModel.FDR <- conf.mat$FDR
                   res$model <- conf.mat$model
                   res$modelConfMatrix <- table(classes, predClasses)
                   res$initModel <- classifier1[1]
                   res$classifier <- classifier2[1]
               }
           }

       # -------------------- To search for a cutpoint --------------------- #
           if (find.cut) {
               if (is.null(classifier2)) classifier2 <- classifier1
               conf.mat <- evaluateDIMPclass(LR, column = column,
                                         control.names = "ctrl",
                                         treatment.names = "treat",
                                         classifier=classifier1[1], prop=prop,
                                         output = "conf.mat",
                                         num.cores=num.cores,
                                         tolerance=tol,
                                         tasks=tasks, verbose = FALSE, ...)

               post <- predict(object = conf.mat$model, newdata = LR,
                               type = "posteriors")

               cuts <- seq(cut.interval[1], cut.interval[2], cut.incr)
               k = 1; opt <- FALSE; overcut <- FALSE;
               while (k < length(cuts) && !opt && !overcut) {
                   cutpoint <- cutFun(divs = unlist(LR), cuts[k])
                   dmps <- selectDIMP(LR, div.col = div.col,
                                   cutpoint = cutpoint,
                                   tv.col=tv.col, tv.cut=tv.cut)
                   if (length(dmps$ctrl) > 0) {
                       conf.mat <- evaluateDIMPclass(dmps, column = column,
                                               control.names = control.names,
                                               treatment.names=treatment.names,
                                               classifier=classifier2[1],
                                               prop=prop, output = "conf.mat",
                                               num.cores=num.cores,
                                               tasks=tasks,
                                               verbose = FALSE, ...)
                       if (stat == 0) {
                           st <- conf.mat$Performance$overall[1]
                           if (st == 1) opt <- TRUE
                           k <- k + 1
                       }
                       if (is.element(stat, 1:11)) {
                           st <- conf.mat$Performance$byClass[stat]
                           if (st == 1) opt <- TRUE
                           k <- k + 1
                       }
                       if (stat == 12) {
                           st <- conf.mat$FDR
                           if (st == 0) opt <- TRUE
                           k <- k + 1
                       }
                   } else {
                       overcut <- TRUE
                       if (k == 1) {
                           st <- conf.mat$Performance$byClass[stat]
                           cutpoint <- cutFun(divs = unlist(LR), cuts[k])
                       }
                       else cutpoint <- cutFun(divs = unlist(LR), cuts[k - 1])
                   }
               }

               predClasses <- predict(object = conf.mat$model, newdata = LR)

               STAT <- c("Accuracy", "Sensitivity", "Specificity",
                       "Pos Pred Value", "Neg Pred Value","Precision", "Recall",
                       "F1", "Prevalence", "Detection Rate",
                       "Detection Prevalence", "Balanced Accuracy", "FDR")

               res$cutpoint <- cutpoint
               res$testSetPerformance <- conf.mat$Performance
               res$testSetModel.FDR <- conf.mat$FDR
               res$model <- conf.mat$model
               res$modelConfMatrix <- table(classes, predClasses)
               res$initModel <- classifier1[1]
               res$postProbCut <- cuts[k]
               res$classifier <- classifier2[1]
               res$statistic <- STAT[stat + 1]
               res$optStatVal <- st
           }
       }
   }
   return(res)
}
