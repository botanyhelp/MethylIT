#' @rdname divTest
#' @title Group Comparisons of Information Divergences Based on Generalized
#'     Linear Model
#' @description Generalized Linear Model for group comparison of information
#'     divergence variables yielded by MethylIT output. Basically, this a
#'     wrapping function to perform the fitting of generalized linear models
#'     with \code{\link[stats]{glm}} from 'stats' package to any variable of
#'     interest given in GRanges objects of MethylIT output.
#' @details The default parameter setting glm.family = Gamma(link = "log") is
#'     thought to perform the group comparison of the sums of absolute
#'     differences of methylation levels (total variation distance (TVD) at
#'     gene-body DIMPs on DMGs). The sums of Hellinger divergence (HD, at
#'     gene-body DIMPs on DMGs) can be tested with this setting as well. Both
#'     TVD and HD follow asymptotic Chi-square distribution and, consequently,
#'     so do the sum of TVD and the sum of HD.  The Chi-square distribution is
#'     a particular case of Gamma distribution: \cr
#'         \deqn{f(x|a,s) = 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)}
#'     Chi-square density is derived after replacing a = n/2 and s = 2: \cr
#'         \deqn{f(x|n) = 1/(2^(n/2) Γ(n/2)) x^(n/2-1) e^(-x/2)}
#' @param GR GRanges objects including control and treatment samples containing
#'     an information divergence of methylation levels. The names for each
#'     column must coincide with the names given for parameters:
#'     'control.names' and 'treatment.names'.
#' @param control.names Names/IDs of the control samples, which must be
#'     included in the variable GR in a metacolumn.
#' @param treatment.names Names/IDs of the treatment samples, which must be
#'     included in the variable GR in a metacolumn.
#' @param glm.family,link Parameter to be passed to function
#'     \code{\link[stats]{glm}}. A description of the error distribution and
#'     link function to be used in the model. For \code{\link[stats]{glm}} this
#'     can be a character string naming a family function, or the result of a
#'     call to a family function. For \code{\link[stats]{glm}}.fit only the
#'     third option is supported. (See\code{\link[stats]{family}} function).
#'     Default: glm.family=Gamma(link ="log").
#' @param var.weights Logical (default: FALSE). Whether to use group variances
#'     as weights.
#' @param weights An optional list of two numeric vectors of ‘prior weights’ to
#'     be used in the fitting process. One vector of weights for the control
#'     and one for the treatment. Each vector with length equal to length(GR)
#'     (default: NULL). Non-NULL weights can be used to indicate that different
#'     observations have different dispersions (with the values in weights
#'     being inversely proportional to the dispersions).
#' @param varFilter Numeric (default: 0). GLM will be performed only for those
#'     rows (ranges denoting genomic regions) where the group variance is
#'     greater the number specified by varFilter.
#' @param meanFilter Numeric (default: 0). GLM will be performed only for those
#'     rows (ranges denoting genomic regions) where the absolute difference of
#'     group means is  greater the number specified by meanFilter.
#' @param FilterLog2FC if TRUE, the results are filtered using the minimun
#'     absolute value of log2FoldChanges observed to accept that a gene in the
#'     treatment is differentially expressed in respect to the control.
#' @param Minlog2FC minimum logarithm base 2 of fold changes
#' @param divPerBp At least for one group the mean divergence per bp must be
#'     equal to or greater than 'divPerBp' (default divPerBp = 0.001).
#' @param minInd Integer (Default: 2). At least one group must have 'minInd'
#'     individuals with a divergence value greater than zero.
#' @param pAdjustMethod Method used to adjust the results; default: "NULL"
#'     (see \code{\link[stats]{p.adjust}}.methods). The p-value adjustment is
#'     performed using function \code{\link[stats]{p.adjust}}.
#' @param scaling integer (default 1). Scaling factor estimate the
#'     signal density as: scaling * "DIMP-Count-Per-Bp". For example,
#'     if scaling = 1000, then signal density denotes the number of DIMPs in
#'      1000 bp.
#' @param pvalCutOff cutoff used then a p-value adjustment is performed
#' @param saveAll if TRUE all the temporal results that passed filters
#'     'varFilter' and are 'meanFilter' returned. If FALSE, only the
#'     comparisons that passed filters 'varFilter', 'meanFilter', and
#'     pvalue < pvalCutOff or adj.pvalue < pvalCutOff (if pAdjustMethod is not
#'     NULL) are returned.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see
#'     \code{\link[BiocParallel]{bplapply}} function from BiocParallel).
#' @param tasks integer(1). The number of tasks per job.  Value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call to
#'     a function, such as bplapply, bpmapply etc. A task is the division of the
#'     X argument into chunks. When tasks == 0 (default), X is divided as evenly
#'     as possible over the number of workers (see
#'     \code{\link[BiocParallel]{MulticoreParam-class}} from BiocParallel
#'     package).
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... Additional parameters passed to \code{\link[stats]{glm}} function.
#' @importFrom stats glm Gamma anova
#' @importFrom genefilter rowVars
#' @importFrom BiocParallel MulticoreParam bplapply
#' @return The original GRanges object with the columns "beta", "log2FC",
#'     "pvalue", "adj.pval" (if pAdjustMethod requested), "CT.divPerBp" and
#'     "TT.divPerBp" (divergence per base pairs), and "divPerBpVariation added.
#' @examples
#' ## Gene annotation
#' genes <- GRanges(seqnames = "1",
#'                  ranges = IRanges(start = c(3631, 6788, 11649),
#'                                   end = c(5899, 9130, 13714)),
#'                  strand = c("+", "-", "-"))
#' mcols(genes) <- data.frame(gene_id = c("AT1G01010", "AT1G01020",
#'                                        "AT1G01030"))
#' # === The number of cytosine sites to generate ===
#' sites = 11001
#' # == Set a seed for pseudo-random number generation ===
#' set.seed(123)
#' alpha.ct <- 0.09
#' alpha.tt <- 0.2
#' # === Simulate samples ===
#' ref = simulateCounts(num.samples = 2, sites = sites, alpha = alpha.ct,
#'                    beta = 0.5, size = 50, theta = 4.5, sample.ids = "R1")
#'
#' # Control group
#' ctrl = simulateCounts(num.samples = 2, sites = sites, alpha = alpha.ct,
#'                        beta = 0.5, size = 50, theta = 4.5,
#'                        sample.ids = c("C1", "C2"))
#' # Treatment group
#' treat = simulateCounts(num.samples = 2, sites = sites, alpha = alpha.tt,
#'                         beta = 0.5, size = 50, theta = 4.5,
#'                         sample.ids = c("T1", "T2"))
#'
#' #  === Estime Divergences ===
#' HD = estimateDivergence(ref = ref$R1, indiv = c(ctrl, treat),
#'                         Bayesian = TRUE, num.cores = 1L, percentile = 1,
#'                         verbose = FALSE)
#'
#' nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)
#'
#' ## Next, the potential signal can be estimated
#' PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 4, alpha = 0.05)
#'
#' ## The cutpoint estimation used to discriminate the signal from the noise
#' cutpoints <- estimateCutPoint(PS, control.names = c("C1", "C2"),
#'                               treatment.names = c("T1", "T2"),
#'                               div.col = 4, verbose = FALSE)
#' ## DIMPs are selected using the cupoints
#' DIMPs <- selectDIMP(PS, div.col = 9, cutpoint = min(cutpoints$cutpoint))
#'
#' ## Finally DIMPs statistics genes
#' tv_DIMPs = getGRegionsStat(GR = DIMPs, grfeatures = genes, stat = "sum",
#'                            absolute = TRUE, column = 7L)
#'
#' GR_tv_DIMP = uniqueGRanges(tv_DIMPs, type = "equal", chromosomes = "1")
#' colnames(mcols(GR_tv_DIMP)) <-  c("C1", "C2", "T1", "T2")
#'
#' res <- divTest(GR=GR_tv_DIMP, control.names =  c("C1", "C2"),
#'                treatment.names = c("T1", "T2"))
#' @export
divTest <- function(GR, control.names, treatment.names,
                   glm.family=Gamma(link="log"), var.weights = FALSE,
                   weights=NULL, varFilter=0, meanFilter=0, FilterLog2FC=TRUE,
                   Minlog2FC=1, divPerBp=0.001, minInd=2, pAdjustMethod=NULL,
                   scaling=1L, pvalCutOff=0.05, saveAll=FALSE, num.cores=1,
                   tasks=0L, verbose=TRUE, ...) {
   if (!inherits(GR, "GRanges")) {
       stop("* 'GR' must be an object from class 'GRanges'")
   }
   GR=GR[, c(control.names, treatment.names)]
   CT=as.matrix(mcols(GR[, control.names]))
   TT=as.matrix(mcols(GR[, treatment.names]))
   m=as.matrix(mcols(GR))
   # === Pre-filtering data === #
   g1=length(control.names); g2=length(treatment.names)
   # At least one group must have 'minInd' individuals with a divergence value
   # greater than zero.
   idx=which((rowSums(CT > 0) >= max(g1 - 1, minInd)) |
               (rowSums(TT > 0) >= max(g2 - 1, minInd)))
   m=m[idx,]; GR=GR[idx]; CT=CT[idx,]; TT=TT[idx,]
   # filtering based on group variance
   vars=rowVars(m)
   idx=which(vars > varFilter)
   m=m[idx,]; GR=GR[idx]; CT=CT[idx,]; TT=TT[idx,]
   ct.mean=rowMeans(CT)
   tt.mean=rowMeans(TT)
   mean.diff=abs(tt.mean - ct.mean)
   idx=which(mean.diff > meanFilter)
   m=m[idx,]; GR=GR[idx]; CT=CT[idx,]; TT=TT[idx,]

   if (!is.null(divPerBp)) {
     ## For each group the divergence per bp must be equal or greater
     ## than divPerBp
     size <- width(GR)
     CT.divPerBp=scaling * rowMeans(CT)/size
     TT.divPerBp=scaling * rowMeans(TT)/size
     idx <- ((CT.divPerBp >= divPerBp) | (TT.divPerBp >= divPerBp))
     m=m[idx,]; GR=GR[idx]; CT=CT[idx,]; TT=TT[idx,]
     CT.divPerBp=CT.divPerBp[idx]; TT.divPerBp=TT.divPerBp[idx]
   } else stop("*** The divergence signal must be divPerBp >= 0")

   if (verbose) {
       message("*** Number of genes after filtering divergences ", nrow(m))
   }
   # === weights === #
   wg=NULL
   if(var.weights==TRUE && is.null(weights)) {
       ct.var=rowVars(CT)
       tt.var=rowVars(TT)
       wg=cbind(1/(ct.mean + ct.var * ct.mean^2),
                1/(tt.mean + tt.var * tt.mean^2))
       rm(CT, TT, ct.var, tt.var); gc()
   }
   if(!is.null(weights)) {
     ct.var=weights[[1]]
     tt.var=weights[[2]]
     wg=cbind(ct.var, tt.var)
   }

   div = t(m); rm(m, vars, ct.mean, tt.mean, idx, mean.diff); gc()
   group <- factor(c(rep(1, length(control.names)),
                   rep(2, length(treatment.names))))
   group <- relevel(group, ref = "1")

   # ===== Auxiliar function to perform beta regression ===== #
   divglm <- function(x, grp, family, weights, ...) {
       model = tryCatch(glm(x ~ grp, family=family, weights=weights, ...),
               error=function(e) {
                   n=length(grp)
                   x=(x * (n - 1) + 0.5)/n
                   try(glm(x ~ grp, family=family, weights=weights, ...),
                       silent=TRUE)})
       if(!is.null(weights)) {
         n=as.vector(table(grp))
         w=c(rep(weights[1], n[1]), rep(weights[2], n[2]))
         rw = tryCatch(glm(x ~ grp, family=family, weights=w, ...),
                       error=function(e) {
                       n=length(grp)
                       x=(x * (n - 1) + 0.5)/n
                       try(glm(x ~ grp, family=family, weights=w, ...),
                           silent=TRUE)})
         if (!inherits(model, "try-error") && !inherits(rw, "try-error")) {
           midx=which.min(c(AIC(model), AIC(rw)))
           model=list(model, rw)[[midx]]
         }
         if (inherits(model, "try-error") && !inherits(rw, "try-error")) {
           model=rw
         }
       }

       if (!inherits(model, "try-error")) {
           beta=unname(coef(model)[2])
           log2FC=beta/log(2)
           pvalue=anova(model, test =  "Chisq")[2, 5]
       } else {
           beta=NA; log2FC=NA; pvalue=NA
       }
           return(c(beta=beta, log2FC=log2FC, pvalue=pvalue))
   }

   if (verbose)
       cat("* Performing generalized linear regression analysis ... \n")
   nc=ncol(div)
   if (num.cores > 1) {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       if (is.null(wg)) {
           res <- bplapply(1:nc, function(k) {
               divglm(x=div[,k], grp=group, family=glm.family, weights=NULL,
                      ...)},
               BPPARAM=bpparam)
       } else {
           res <- bplapply(1:nc, function(k) {
               divglm(x=div[,k], grp=group, family=glm.family,
                      weights=wg[k,], ...)}, BPPARAM=bpparam)
       }
       res=do.call(rbind, res)
   } else {
       res=matrix(NA, nrow=nc, ncol=3)
       if (is.null(wg)) {
           for (k in 1:nc) {
               res[k,]=divglm(x=div[,k], grp=group, family=glm.family,
                              weights=NULL)}
       } else {
           for (k in 1:nc) {
               res[k,]=divglm(x=div[,k], grp=group, family=glm.family,
                       weights=wg[k,])
           }
       }
       colnames(res) <- c("beta", "log2FC", "pvalue")
   }
   res = data.frame(res)
   res$CT.divPerBp <- CT.divPerBp
   res$TT.divPerBp <- TT.divPerBp
   res$divPerBpVariation <- (TT.divPerBp - CT.divPerBp)

   # === Filtering results ===
   if (FilterLog2FC) fidx <- which(abs(res$log2FC) > Minlog2FC)
   if (!is.null(pAdjustMethod)) {
       res$adj.pval < rep(NA, length(res$pvalue))
       # pvalue adjustment only for the selected comparisons: "fidx"
       if (FilterLog2FC) {
           res$adj.pval[fidx] <- p.adjust(res$pvalue[fidx],
                                  method=pAdjustMethod)
       } else
         # pvalue adjustment for the comparisons
         res$adj.pval <- p.adjust(res$pvalue, method=pAdjustMethod)
   }
   # Only results with "adj.pval < pvalCutOff" are reported if saveAll=FALSE
   if (!is.null(pvalCutOff) && !saveAll) {
      if (!is.null(pAdjustMethod)) {
           idx=which(res$adj.pval < pvalCutOff)
           res <- res[idx, ]
           GR=GR[idx]
      } else {
           idx=which(res$pvalue < pvalCutOff)
           res <- res[idx, ]
           GR=GR[idx]
      }
   }
   # Only results with "log2FC > Minlog2FC" are reported if saveAll=FALSE
   if (FilterLog2FC && !saveAll) {
       idx=which(abs(res$log2FC) > Minlog2FC)
       res <- res[idx, ]
       GR=GR[idx]
   }

  mcols(GR) <- data.frame(GR@elementMetadata, res)
  return(sortBySeqnameAndStart(GR))
}
