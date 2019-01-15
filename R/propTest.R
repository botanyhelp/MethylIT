#' @rdname propTest
#' @title Beta Regression for methylation levels and rates
#' @description Beta Regression analysis for treatment versus control group
#'     comparison of methylation levels, appends three new metacolumns "beta",
#'     "log2FC", "pvalue" to the provided GRanges argument
#' @details Beta Regression analysis for group comparison of methylation levels
#'     is performed using the function \code{\link[betareg]{betareg}}.
#' @param GR GRanges objects including control and treatment samples containing
#'     the methylation levels. The name for each column must coincide with the
#'     names given for parameters: 'control.names' and 'treatment.names'.
#' @param control.names Names/IDs of the control samples, which must be include
#'     in the variable GR at the metacolumn.
#' @param treatment.names Names/IDs of the treatment samples, which must be
#'     included in the variable GR at the metacolumn.
#' @param link Parameter to be passed to function 'betareg' from package
#'     'betareg'. character specification of the link function in the mean model
#'     (mu). Currently, "logit", "probit", "cloglog", "cauchit", "log",
#'     "loglog" are supported. Alternatively, an object of class "link-glm"
#'     can be supplied (see \code{\link[betareg]{betareg}}).
#' @param type Parameter to be passed to function 'betareg' from package
#'     'betareg'. A character specification of the type of estimator. Currently,
#'     maximum likelihood ("ML"), ML with bias correction ("BC"), and ML with
#'     bias reduction ("BR") are supported.
#' @param tv.cut A cutoff for the total variation distance (TVD; absolute value
#'     of methylation levels differences) estimated at each site/range as the
#'     difference of the group means of methylation levels. If tv.cut is
#'     provided, then sites/ranges k with abs(TV_k) < tv.cut are removed before
#'     to perform the regression analysis. Its value must be NULL or a number
#'     0 < tv.cut < 1.
#' @param indvPerGrp An integer number giving the minimum number of individuals
#'     per group at each site/region. Default 2.
#' @param FilterLog2FC if TRUE, the results are filtered using the minimun
#'     absolute value of log2FoldChanges observed to accept that a gene in the
#'     treatment is differentially expressed in respect to the control.
#' @param pAdjustMethod method used to adjust the results; default: BH
#' @param pvalCutOff cutoff used, then a p-value adjustment is performed. If
#'     NULL all the reported p-values are for testing.
#' @param Minlog2FC minimum logarithm base 2 of fold changes.
#' @param saveAll if TRUE all the temporal results are returned.
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
#' @importFrom betareg betareg
#' @importFrom BiocParallel MulticoreParam bplapply
#' @return The original GRanges object with the columns "beta", "log2FC",
#'     "pvalue", and TV added.
#' @examples
#' num.cyt <- 11001 # Number of cytosine position with methylation call
#' max.cyt = 14000
#' ## Gene annotation
#' genes <- GRanges(seqnames = "1",
#'                  ranges = IRanges(start = c(3631, 6788, 11649),
#'                                   end = c(5899, 9130, 13714)),
#'                  strand = c("+", "-", "-"))
#' mcols(genes) <- data.frame(gene_id = c("AT1G01010", "AT1G01020",
#'                                        "AT1G01030"))
#'
#' set.seed(123) #'#' To set a seed for random number generation
#' ## GRanges object of the reference with methylation levels in
#' ## its meta-column
#' Ref <- makeGRangesFromDataFrame(
#'   data.frame(chr = '1',
#'              start = 3000:max.cyt,
#'              end = 3000:max.cyt,
#'              strand = '*',
#'              p1 = rbeta(num.cyt, shape1 = 1, shape2 = 1.5)),
#'   keep.extra.columns = TRUE)
#'
#' ## List of Granges objects of individuals methylation levels
#' Indiv <- GRangesList(
#'   sample11 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 3000:max.cyt,
#'                end = 3000:max.cyt,
#'                strand = '*',
#'                p2 = rbeta(num.cyt, shape1 = 1.5, shape2 = 2)),
#'     keep.extra.columns = TRUE),
#'   sample12 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 3000:max.cyt,
#'                end = 3000:max.cyt,
#'                strand = '*',
#'                p2 = rbeta(num.cyt, shape1 = 1.6, shape2 = 2.1)),
#'     keep.extra.columns = TRUE),
#'   sample21 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 3000:max.cyt,
#'                end = 3000:max.cyt,
#'                strand = '*',
#'                p2 = rbeta(num.cyt, shape1 = 10, shape2 = 4)),
#'     keep.extra.columns = TRUE),
#'   sample22 = makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                start = 3000:max.cyt,
#'                end = 3000:max.cyt,
#'                strand = '*',
#'                p2 = rbeta(num.cyt, shape1 = 11, shape2 = 4)),
#'     keep.extra.columns = TRUE))
#' ## To estimate Hellinger divergence using only the methylation levels.
#' HD <- estimateDivergence(ref = Ref, indiv = Indiv, meth.level = TRUE,
#'                          columns = 1)
#' ## To perform the nonlinear regression analysis
#' nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)
#'
#' ## Next, the potential signal can be estimated
#' PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 4, alpha = 0.05)
#'
#' ## The cutpoint estimation used to discriminate the signal from the noise
#' cutpoints <- estimateCutPoint(PS, control.names = c("sample11", "sample12"),
#'                               treatment.names = c("sample21", "sample22"),
#'                               div.col = 4, verbose = TRUE)
#' ## DIMPs are selected using the cupoints
#' DIMPs <- selectDIMP(PS, div.col = 4, cutpoint = min(cutpoints$cutpoint))
#'
#' ## Finally DIMPs statistics genes
#' p_DIMPs = getGRegionsStat(GR = DIMPs, grfeatures = genes, stat = "mean",
#'                           prob = TRUE, column = 2L)
#'
#' GR_p_DIMP = uniqueGRanges(p_DIMPs, type = "equal", chromosomes = "1")
#' colnames(mcols(GR_p_DIMP)) <-  c("sample11", "sample12", "sample21",
#'                                 "sample22")
#' names(GR_p_DIMP) <- genes$gene_id
#'
#' ## Group differences between methylation levels
#' propTest(GR = GR_p_DIMP, control.names = c("sample11", "sample12"),
#'          treatment.names = c("sample21", "sample22"))
#' @export
propTest <- function(GR, control.names, treatment.names, link="logit",
                     type="ML", tv.cut=NULL, indvPerGrp = 0, FilterLog2FC=TRUE,
                     pAdjustMethod="BH", pvalCutOff=0.05, Minlog2FC=0.5,
                     saveAll=FALSE, num.cores=1, tasks=0L, verbose=TRUE) {
   if (!inherits(GR, "GRanges")) {
    stop("* 'GR' must be an object from class 'GRanges'")
   }
   GR=GR[, c(control.names, treatment.names)]
   prop=as.matrix(mcols(GR))
   g1=match(control.names, c(control.names, treatment.names))
   g2=match(treatment.names, c(control.names, treatment.names))
   idx <- (rowSums(prop[, g1] > 0) >= indvPerGrp |
           rowSums(prop[, g2] > 0) >= indvPerGrp)
   GR = GR[idx]
   prop <- prop[idx, ]
   GR$TV <- rowMeans(prop[, g2]) - rowMeans(prop[, g1])
   if (!is.null(tv.cut)) {
     idx.tv <- (abs(GR$TV) >= tv.cut)
     if(sum(idx.tv) > 0) {
       idx.tv = which(idx.tv)
       prop <- prop[idx.tv, ]
     } else stop("*** tv.cut is higher than data TV values")
   }

   prop <- t(prop)
   group <- factor(c(rep(1, length(g1)),
                   rep(2, length(g2))))
   group <- relevel(group, ref = "1")

   # Auxiliar function to perform beta regression
   betaReg <- function(x, grp, link, type) {
       r = tryCatch(betareg(x ~ grp, link=link, type=type),
               error=function(e) {
               n=length(grp)
               x=(x * (n - 1) + 0.5)/n
               try(betareg(x ~ grp, link=link, type=type),
                   silent=TRUE)})
   if (inherits(r, "try-error")) {
       r=try(glm(x ~ grp, family=quasibinomial(link=link)),
           silent=TRUE)
   }

   if (!inherits(r, "try-error")) {
       beta=unname(coef(r)[2])
       log2FC=beta/log(2)
       if (inherits(r, "betareg")) {
           pvalue=unname(coefficients(summary(r))$mean[2, 4])
       } else pvalue=coef(summary(r))[2,4]
   } else {
       beta=NA; log2FC=NA; pvalue=NA
   }
   return(c(beta=beta, log2FC=log2FC, pvalue=pvalue))
   }

   if (verbose) cat("*** Performing beta-regression analysis ... \n")
   n=ncol(prop)
   if(num.cores > 1) {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       res <- bplapply(1:n, function(k) {
           betaReg(x=prop[,k], grp=group, link=link, type=type)},
           BPPARAM=bpparam)
   res=do.call(rbind, res)
   res = data.frame(res)
   } else {
       res=matrix(NA, nrow=n, ncol=3)
       for (k in 1:n) {
           res[k,]=betaReg(x=prop[,k], grp=group, link=link, type=type)
       }
       res = data.frame(res)
       colnames(res) <- c("beta", "log2FC", "pvalue")
   }

   if (all(is.na(res))) stop("*** The regression model does not fit your data")

   if (saveAll && !is.null(tv.cut)) {
       l = length(GR)
       GR$beta <- rep(NA, l)
       GR$beta[idx.tv] <- res$beta
       GR$log2FC <- rep(NA, l)
       GR$log2FC[idx.tv] <- res$log2FC
       GR$pvalue <- rep(NA, l)
       GR$pvalue[idx.tv] <- res$pvalue
   } else {
       if (!is.null(tv.cut)) {
           GR = GR[idx.tv]
           mcols(GR) <- data.frame(GR@elementMetadata, res)
       }
       if (is.null(tv.cut)) mcols(GR) <- data.frame(GR@elementMetadata, res)
   }

   if (!saveAll) GR <- GR[!is.na(GR$log2FC)]
   if (FilterLog2FC && is.null(pvalCutOff) && !saveAll) {
       GR <- GR[abs(GR$log2FC) > Minlog2FC]
   }
   GR$adj.pval <- rep(NA, length(GR))
   if (!is.null(pvalCutOff) && !saveAll) {
       if (FilterLog2FC) GR <- GR[abs(GR$log2FC) > Minlog2FC, ]
       GR$adj.pval <- p.adjust(GR$pvalue, method=pAdjustMethod)
       GR <- GR[GR$adj.pval < pvalCutOff, ]
   } else {
     if (!is.null(pvalCutOff) && saveAll && !FilterLog2FC) {
       idx = which(!is.na(GR$adj.pval))
       GR$adj.pval[idx] <- p.adjust(GR$pvalue[idx], method=pAdjustMethod)
     }
     if (FilterLog2FC) {
       idx <- which(abs(GR$log2FC) > Minlog2FC)
       pval <- GR$pvalue[idx]
       GR$adj.pval[idx] <- p.adjust(pval, method=pAdjustMethod)
     } else {
       GR$adj.pval <- p.adjust(GR$pvalue, method=pAdjustMethod)
     }
   }
   return(GR)
}
