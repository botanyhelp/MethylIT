#' @rdname estimateCutPoint
#'
#' @title Estimate cutpoints to distinguish the treatment methylation signal
#'     from the control
#' @description Given a list of two GRanges objects, control and treatment,
#'     carrying the potential signals (prior classification) from controls and
#'     treatments in terms of an information divergence
#'     (given the meta-columns), the function estimates the cutpoints of the
#'     control group versus treatment group.
#' @details The function performs an estimation of the optimal cutpoint for the
#'     classification of the differentially methylated (cytosines) positions
#'     into two classes: DMPs from control and DMPs from treatment. The simplest
#'     approach to estimate the cutpoint  is based on the application of Youden
#'     Index. More complexes approach based in several machine learning model
#'     are provided as well.
#'
#'     Results of the classification perfomance resulting from the estimated
#'     cutpoint are normally given, with the exception of those extreme
#'     situations where the statistics to evaluate performance cannot be
#'     estimated. More than one classifier model can be applied. For example,
#'     one classifier (logistic model) can be used to esitmate the posterior
#'     classification probabilities of DMP into those from control and those
#'     from treatment. This probabilities are then used to estimate the cutpoint
#'     in range of values from, say, 0.5 to 0.8. Next, a different classifier
#'     can be used to evaluate the classification performance. Different
#'     classifier models would yield different performances. Models are returned
#'     and can be used in futher prediction with new datasets from the same
#'     batch experiment. This is a machine learnig approach to discriminate the
#'     biological regulatory signal naturally generated in the control from that
#'     one induced by the treatment.
#'
#' @param LR An object from 'pDMP' class. This obejct is previously obtained
#'     with function \code{\link{getPotentialDIMP}}.
#' @param control.names,treatment.names Names/IDs of the control and
#'     treatment samples, which must be include in the variable LR.
#' @param simple Logic (default, TRUE). If TRUE, then Youden Index is used to
#'     estimate the cutpoint.
#' @param column a logical vector for column names for the predictor variables
#'     to be used: Hellinger divergence "hdiv", total variation "TV",
#'     probability of potential DIMP "wprob", and the relative cytosine site
#'     position "pos" in respect to the chromosome where it is located. The
#'     relative position is estimated as (x - x.min)/(x.max - x), where x.min
#'     and x.max are the maximum and minimum for the corresponding chromosome,
#'     repectively. If "wprob = TRUE", then Logarithm base-10 of "wprob" will
#'     be used as predictor in place of "wprob".
#' @param classifier1,classifier2 Classification model to use. Option "logistic"
#'     applies a logistic regression model; option "lda" applies a Linear
#'     Discriminant Analysis (LDA); "qda" applies a Quadratic Discriminant
#'     Analysis (QDA), "pca.logistic" applies logistic regression model using
#'     the Principal Component (PCs) estimated with Principal Component Analysis
#'     (PCA) as predictor variables. pca.lda" applies LDA using PCs as predictor
#'     variables, and the option "pca.qda" applies a Quadratic Discriminant
#'     Analysis (QDA) using PCs as predictor variables. If classifier2 is not
#'     NULL, then it will be used to evaluate the classification performance,
#'     and the corresponding best fitted model will be returned.
#' @param tv.cut A cutoff for the total variation distance to be applied to each
#'     site/range. Only sites/ranges \emph{k} with \eqn{TVD_{k} > tv.cut} are
#'     are used in the analysis. Its value must be a number.
#'     \eqn{0 < tv.cut < 1}. Default is \eqn{tv.cut = 0.25}.
#' @param div.col Column number for divergence variable for which the estimation
#'     of the cutpoint will be performed.
#' @param clas.perf Logic. Whether to evaluate the classification performance
#'     for the estimated cutpoint using a model classifier when 'simple=TRUE'.
#'     Default, FALSE.
#' @param post.cut If 'simple=FALSE', this is posterior probability to dicide
#'     whether a DMPs belong to treatment group. Default \emph{post.cut} = 0.5.
#' @param prop Proportion to split the dataset used in the logistic regression:
#'     group versus divergence (at DIMPs) into two subsets, training and
#'     testing.
#' @param n.pc Number of principal components (PCs) to use if the classifier is
#'     not 'logistic'. In the current case, the maximun number of PCs is 4.
#' @param find.cut Logic. Wether to search for an optimal cutoff value to
#'     classify DMPs based on given specifications.
#' @param cut.interval 0 < \emph{cut.interval} < 0.1. If \emph{find.cut} = TRUE, the
#'     interval of treatment group posterior probabilities where to search for a
#'     cutpoint. Deafult \emph{cut.interval} = c(0.5, 0.8).
#' @param cut.incr 0 < \emph{cut.incr} < 0.1. If \emph{find.cut} = TRUE, the sucesive
#'     increamental values runing on the interval \emph{cut.interval}. Deafult,
#'     \emph{cut.incr} = 0.01.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param stat An integer number indicating the statistic to be used in the
#'     testing when \emph{find.cut} = TRUE. The mapping for statistic names are:
#'     \itemize{
#'         \item 0 = "Accuracy"
#'         \item 1 = "Sensitivity"
#'         \item 2 = "Specificity"
#'         \item 3 = "Pos Pred Value"
#'         \item 4 = "Neg Pred Value"
#'         \item 5 = "Precision"
#'         \item 6 = "Recall"
#'         \item 7 = "F1"
#'         \item 8 = "Prevalence"
#'         \item 9 = "Detection Rate"
#'         \item 10 = "Detection Prevalence"
#'         \item 11 = "Balanced Accuracy"
#'         \item 12 = FDR
#'    }
#' @return Depending the parameter setting will return the following list with
#'     elements:
#'     \enumerate{
#'         \item cutpoint: Cutpoint estimated.
#'         \item testSetPerformance: Performance evaluation on the test set.
#'         \item testSetModel.FDR: False discovery rate on the test set.
#'         \item model: Model used in the performance evaluation.
#'         \item modelConfMatrix: Confusion matrix for the whole dataset derived
#'               applying the model classifier used in the performance
#'               evaluation.
#'         \item initModel: Initial classifier model applied to estimate
#'               posterior classifications used in the cutpoint estimation.
#'         \item postProbCut: Posterior probability used to estimate the
#'               cutpoint
#'         \item classifier: Name of the model classifier used in the
#'               performance evaluation.
#'         \item statistic: Name of the performance statistic used to find the
#'               cutpoint when find.cut =  TRUE.
#'         \item optStatVal: Value of the performance statistic at the cutpoint.
#'     }
#'
#' @examples
#' set.seed(123) ## To set a seed for random number generation
#' ## GRanges object of the reference with methylation levels in
#' ## its matacolumn
#' num.points <- 5000
#' Ref <- makeGRangesFromDataFrame(
#'   data.frame(chr = '1',
#'              start = 1:num.points,
#'              end = 1:num.points,
#'              strand = '*',
#'              p1 = rbeta(num.points, shape1 = 1, shape2 = 1.5)),
#'   keep.extra.columns = TRUE)
#'
#' ## List of Granges objects of individuals methylation levels
#' Indiv <- GRangesList(
#'   sample11 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 1:num.points,
#'                end = 1:num.points,
#'                strand = '*',
#'                p2 = rbeta(num.points, shape1 = 1.5, shape2 = 2)),
#'     keep.extra.columns = TRUE),
#'   sample12 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 1:num.points,
#'                end = 1:num.points,
#'                strand = '*',
#'                p2 = rbeta(num.points, shape1 = 1.6, shape2 = 2)),
#'     keep.extra.columns = TRUE),
#'   sample21 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 1:num.points,
#'                end = 1:num.points,
#'                strand = '*',
#'                p2 = rbeta(num.points, shape1 = 40, shape2 = 4)),
#'     keep.extra.columns = TRUE),
#'   sample22 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 1:num.points,
#'                end = 1:num.points,
#'                strand = '*',
#'                p2 = rbeta(num.points, shape1 = 41, shape2 = 4)),
#'     keep.extra.columns = TRUE))
#' ## To estimate Hellinger divergence using only the methylation levels.
#' HD <- estimateDivergence(ref = Ref, indiv = Indiv, meth.level = TRUE,
#'                          columns = 1)
#' ## To perform the nonlinear regression analysisx
#' nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)
#'
#' ## Next, the potential signal can be estimated
#' PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 4, alpha = 0.05)
#'
#' cutpoint <- estimateCutPoint(LR = PS, simple = TRUE, find.cut = FALSE,
#'                              column = c(hdiv = TRUE, TV = TRUE,
#'                                         wprob = TRUE, pos = TRUE),
#'                              interaction = "hdiv:TV", clas.perf = FALSE,
#'                              control.names = c("sample11", "sample12"),
#'                              treatment.names = c("sample21", "sample22"))
#' @importFrom S4Vectors mcols
#' @importFrom caret confusionMatrix
#' @export
estimateCutPoint <- function(LR, control.names, treatment.names, simple = TRUE,
                       column=c(hdiv=TRUE, TV=TRUE, wprob=FALSE, pos=FALSE),
                       classifier1=c("logistic", "pca.logistic", "lda",
                                   "qda","pca.lda", "pca.qda"),
                       classifier2=NULL, tv.cut = 0.25, div.col = NULL,
                       clas.perf = FALSE, post.cut = 0.5, prop=0.6, n.pc=1,
                       find.cut=FALSE, cut.interval = c(0.5, 0.8),
                       cut.incr = 0.01, stat = 1, num.cores=1L, tasks=0L, ...) {
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
   # -----------------------Divergences  are positives ------------------------#
   # In case that TV column would be used as source to get TVD
   vn <- c("hdiv", "TV")
   vn <- vn[match(TRUE, column[vn])]
   if (is.null(div.col)) div.col <- vn

   sn <- names(LR)
   if (any(unlist(
       lapply(LR, function(GR) min(mcols(GR[, div.col])[, 1],
                                   na.rm = TRUE))) < 0)) {
        LR <- lapply(LR, function(GR) {
           mcols(GR[, div.col])[, 1] <- abs(mcols(GR[, div.col])[, 1])
           return(GR)
       }, keep.attr = TRUE)
   }

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
       return(as.numeric(names(idx)))
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
       LR <- try(LR[c(control.names, treatment.names)], silent=TRUE)
   if (inherits(LR, "try-error"))
       stop("List's names does not match control & treatment names")

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
       classes <- factor(classes, levels = c("CT", "TT"))

       if (simple) {
           cutpoint <- roc(dt = infDiv(LR, div.col = div.col))
           predClasses <- unlist(LR)$hdiv > cutpoint
           predClasses[ predClasses == TRUE ] <- "TT"
           predClasses[ predClasses == FALSE ] <- "CT"
           predClasses <- factor(predClasses, levels = c("CT", "TT"))
           conf.mat <- confusionMatrix(data=predClasses, reference=classes,
                                       positive="TT")

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
                                           tasks=tasks, verbose = FALSE, ...)

               predClasses <- predict(object = conf.mat$model, newdata = LR)
               predClasses <- factor(predClasses, levels = c("CT", "TT"))
               conf.matrix <- confusionMatrix(data=predClasses,
                                               reference=classes,
                                               positive="TT")

               res$cutpoint <- cutpoint
               res$testSetPerformance <- conf.mat$Performance
               res$testSetModel.FDR <- conf.mat$FDR
               res$model <- conf.mat$model
               res$modelConfMatrix <- conf.matrix
               res$initModel <- "Youden Index"
               res$classifier <- classifier1[1]
           } else {
               res$cutpoint <- cutpoint
               res$modelConfMatrix <- conf.mat
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
                   predClasses <- factor(predClasses, levels = c("CT", "TT"))
                   conf.matrix <- confusionMatrix(data=predClasses,
                                                  reference=classes,
                                                  positive="TT")

                   res$postProbCut <- 0.5
                   res$cutpoint <- cutFun(divs = unlist(LR), 0.5)
                   res$testSetPerformance <- conf.mat$Performance
                   res$testSetModel.FDR <- conf.mat$FDR
                   res$model <- conf.mat$model
                   res$modelConfMatrix <- conf.matrix
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
                                               tasks=tasks, verbose = FALSE,
                                               ...)

                   predClasses <- predict(object = conf.mat$model, newdata = LR)
                   predClasses <- factor(predClasses, levels = c("CT", "TT"))
                   conf.matrix <- confusionMatrix(data=predClasses,
                                                  reference=classes,
                                                  positive="TT")

                   res$cutpoint <- cutpoint
                   res$postProbCut <- post.cut
                   res$testSetPerformance <- conf.mat$Performance
                   res$testSetModel.FDR <- conf.mat$FDR
                   res$model <- conf.mat$model
                   res$modelConfMatrix <- conf.matrix
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
               predClasses <- factor(predClasses, levels = c("CT", "TT"))
               conf.matrix <- confusionMatrix(data=predClasses,
                                              reference=classes,
                                              positive="TT")

               STAT <- c("Accuracy", "Sensitivity", "Specificity",
                       "Pos Pred Value", "Neg Pred Value","Precision", "Recall",
                       "F1", "Prevalence", "Detection Rate",
                       "Detection Prevalence", "Balanced Accuracy", "FDR")

               res$cutpoint <- cutpoint
               res$testSetPerformance <- conf.mat$Performance
               res$testSetModel.FDR <- conf.mat$FDR
               res$model <- conf.mat$model
               res$modelConfMatrix <- conf.matrix
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
