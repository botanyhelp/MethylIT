#' @rdname rmstGR
#'
#' @title Root Mean Square Test for metadata in a list of GRanges objects
#' @description Count data in MethylIT pipeline is carried in GRanges objects.
#'     This function provides a shortcut to apply the parametric Bootstrap
#'     of 2x2 Contingency independence test, which is implemented in function
#'     \code{\link[MethylIT]{bootstrap2x2}}. The  goodness of fit statistic is
#'     the root-mean-square statistic (RMST) or Hellinger divergence, as
#'     proposed by Perkins et al. [1, 2]. Hellinger divergence (HD) is computed
#'     as proposed in [3].
#' @details Samples from each group are pooled according to the statistic
#'     selected (see parameter pooling.stat) and a unique GRanges object is
#'     created with the methylated and unmethylated read counts for each group
#'     (control and treatment) in the metacolumn. So, a contingency table can be
#'     built for range from GRanges object.
#' @param LR A list of GRanges, a GRangesList, a CompressedGRangesList object.
#'     Each GRanges object from the list must have two columns: methylated
#'     (mC) and unmethylated (uC) counts. The name of each element from the
#'     list must coincide with a control or a treatment name.
#' @param control.names Names/IDs of the control samples, which must be include
#'     in the variable GR at the metacolumn.
#' @param treatment.names Names/IDs of the treatment samples, which must be
#'     included in the variable GR at the metacolumn.
#' @param stat Statistic to be used in the testing: 'RMST' or 'HD'
#' @param pooling.stat statistic used to estimate the methylation pool: row sum,
#'     row mean or row median of methylated and unmethylated read counts across
#'     individuals. If the number of control samples is greater than 2 and
#'     pooling.stat is not NULL, then they will pooled. The same for treatment.
#'     Otherwise, all the pairwise comparisons will be done.
#' @param tv.cut A cutoff for the total variation distance (TVD; absolute value
#'     of methylation levels differences) estimated at each site/range as the
#'     difference of the group means of methylation levels. If tv.cut is
#'     provided, then sites/ranges k with abs(TV_k) < tv.cut are removed before
#'     to perform the regression analysis. Its value must be NULL or a number
#'     0 < tv.cut < 1.
#' @param num.permut Number of permutations.
#' @param pAdjustMethod method used to adjust the results; default: BH
#' @param pvalCutOff cutoff used when a p-value adjustment is performed
#' @param saveAll if TRUE all the temporal results are returned
#' @param mc.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see bpapply function from
#'     BiocParallel).
#' @param tasks integer(1). The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... Additional parameters for function
#'     \code{\link[uniqueGRanges]{MethylIT}}.
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @return A GRanges object with the original sample counts, bootstrap p-value
#'     probability, total variation (difference of methylation levels), and
#'     p-value adjusment.
#'
#' @examples
#' #' A list of GRanges
#' set.seed(123)
#' sites = 50
#' data <- list(
#'   C1 = data.frame(chr = "chr1", start = 1:sites,
#'                   end = 1:sites,strand = '*',
#'                   mC = rnbinom(size = 8, mu = 3, n = sites),
#'                   uC = rnbinom(size = 50, mu = 10, n = sites)),
#'   C2 = data.frame(chr = "chr1", start = 1:sites,
#'                   end = 1:sites, strand = '*',
#'                   mC = rnbinom(size = 8, mu = 3, n = sites),
#'                   uC = rnbinom(size = 50, mu = 10, n = sites)),
#'   T1 = data.frame(chr = "chr1", start = 1:sites,
#'                   end = 1:sites,strand = '*',
#'                   mC = rnbinom(size = 50, mu = 10, n = sites),
#'                   uC = rnbinom(size = 10, mu = 10, n = sites)),
#'   T2 = data.frame(chr = "chr1", start = 1:sites,
#'                   end = 1:sites, strand = '*',
#'                   mC = rnbinom(size = 50, mu = 10, n = sites),
#'                   uC = rnbinom(size = 5, mu = 10, n = sites)))
#' #' Transforming the list of data frames into a list of GRanges objects
#' data = lapply(data,
#'               function(x)
#'                 makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))
#'
#' rmstGR(LR = data, control.names = c("C1", "C2"),
#'        treatment.names = c("T1", "T2"),
#'        tv.cut = 0.25, num.permut = 200, pAdjustMethod="BH",
#'        pvalCutOff = 0.05, num.cores = 4L, verbose=TRUE)
#' @references
#'     1. Perkins W, Tygert M, Ward R. and Classical Exact Tests Often
#'         Wildly Misreport Significance; the Remedy Lies in Computers
#'         [Internet]. Uploaded to ArXiv. 2011. Report No.: arXiv:1108.4126v2.
#'     2. Perkins, W., Tygert, M. & Ward, R. Computing the confidence levels
#'         for a root-mean square test of goodness-of-fit. 217, 9072-9084
#'         (2011).
#'     3. Basu, A., Mandal, A. & Pardo, L. Hypothesis testing for two discrete
#'         populations based on the Hellinger distance. Stat. Probab. Lett. 80,
#'         206-214 (2010).
#'
#' @seealso \code{\link[MethylIT]{FisherTest}}
#' @export
rmstGR <- function(LR, control.names, treatment.names, stat="rmst",
                   pooling.stat = "sum", tv.cut=NULL, num.permut=100,
                   pAdjustMethod="BH", pvalCutOff=0.05, saveAll=FALSE,
                   num.cores=1L, tasks=0L, verbose=TRUE, ...) {
   if (inherits(LR, "list")) LR <- try(as(LR, "GRangesList"))
   if (inherits(LR, "try-error"))
       stop("LR is not list of GRanges objects")

   if (!is.element(class(LR), c("GRangesList", "CompressedGRangesList"))) {
       stop("LR is not GRangesList, CompressedGRangesList, GRanges
           or a list of GRanges object")
   }
   LR = try(LR[c(control.names, treatment.names)], silent=TRUE)
   if (inherits(LR, "try-error"))
     stop("List's names does not match control & treatment names")

   # === Auxiliar function to perform RMST ===
   rmst <- function(LR, ctrl.ns, treat.ns) {
       if(length(ctrl.ns) > 2) {
           if (verbose)
               cat("*** # of control samples is > 2. Pooling samples ...", "\n")
               ctrl <- poolFromGRlist(LR[ctrl.ns], stat=pooling.stat,
                               num.cores=num.cores, verbose=verbose)
       } else ctrl <- LR[ctrl.ns][[1]]

       if(length(treat.ns) > 2) {
           if (verbose)
               cat("*** # of treatment samples is > 2. Pooling samples ...",
                   "\n")
       treat <- poolFromGRlist(LR[treat.ns], stat=pooling.stat,
                               num.cores=num.cores, verbose=verbose)
       } else  treat <- LR[treat.ns][[1]]

       GR = uniqueGRanges(list(ctrl, treat), ...)
       count.matrix = as.matrix(mcols(GR))
       p1 <- count.matrix[, 1:2]
       p1 <- p1[, 1]/rowSums(p1)
       p2 <- count.matrix[, 3:4]
       p2 <- p2[, 1]/rowSums(p2)
       TV = p2 - p1; rm(p1, p2); gc()

       if (!is.null(tv.cut)) {
           idx <- which(abs(TV) >= tv.cut)
           count.matrix = as.matrix(mcols(GR[idx]))
       } else count.matrix = as.matrix(mcols(GR))
       sites = nrow(count.matrix)
       GR$TV <- TV
       GR$pvalue <- rep(1, length(GR))
       GR$adj.pval <- rep(1, length(GR))
       count.matrix=split(count.matrix, row(count.matrix))

       if (verbose)
           cat("*** Performing bootstrap analysis for contingence table... \n
               # of sites after filtering: ", sites, "\n")
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
           pvals <- unname(unlist(bplapply(count.matrix, function(v) {
               m=matrix(as.integer(v), 2, byrow = TRUE)
               bootstrap2x2(x=m, stat=stat, num.permut=num.permut)},
               BPPARAM=bpparam)))

       if (!is.null(tv.cut) && !saveAll) {
           GR <- GR[idx]
           GR$pvalue <- pvals
           GR$adj.pval <- p.adjust(pvals, method=pAdjustMethod)
       } else {
           if (!is.null(tv.cut) && saveAll) {
              GR$pvalue[idx] <- pvals
              GR$adj.pval[idx] <- p.adjust(pvals, method=pAdjustMethod)
           }
           if (is.null(tv.cut) && saveAll) {
              GR$pvalue <- pvals
              GR$adj.pval <- p.adjust(pvals, method=pAdjustMethod)
           }
       }

       if (!is.null(pvalCutOff) && !saveAll) {
           GR <- GR[ GR$adj.pval < pvalCutOff ]
       }

       return(GR)
   }

   if (!is.null(pooling.stat) && length(LR) > 2) {
       res <- rmst(LR=LR, ctrl.ns=control.names, treat.ns=treatment.names)
   } else {
       ctrl.ns = control.names
       treat.ns = treatment.names
       res = list()
       i = 1
       test.name = c()
       for(k in 1:length(ctrl.ns)){
           for (j in 1:length(treat.ns)) {
               test.name = c(test.name, paste0(ctrl.ns[k], ".", treat.ns[j]))
               if (verbose)
                   cat("*** Testing",
                       paste0(ctrl.ns[k], " versus ", treat.ns[j]), "\n")
               res[[i]]=rmst(LR=LR, ctrl.ns=control.names[j],
                               treat.ns=treatment.names[k])
               i = i + 1
           }
       }
       names(res) <- test.name
   }
   return(res)
}
