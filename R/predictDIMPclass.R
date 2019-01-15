#' @name predictDIMPclass
#' @rdname predictDIMPclass
#' @title Predict DIMP class
#' @description This function classify each DIMP as a control or a treatment
#'   DIMP
#' @details Predictions only makes sense if the query DIMPs belong to same
#'   methylation context and derive from an experiment accomplished under the
#'   same condition set for the DIMPs used to build the model.
#' @param LR A list of GRanges objects obtained through the through MethylIT
#'   downstream analysis. Basically, this object is a list of GRanges containing
#'   only differentially informative position (DIMPs). The metacolumn of each
#'   GRanges must contain the columna: Hellinger divergence "hdiv", total
#'   variation "TV", the probability of potential DIMP "wprob", which naturally
#'   are added in the downstream analysis of MethylIT.
#' @param model A classifier model obtained with the function
#' 'evaluateDIMPclass'.
#' @param conf.matrix Optional. Logic, whether a confusion matrix should be
#'   returned (default, FALSE, see below).
#' @param control.names Optional. Names/IDs of the control samples, which must
#'   be include in thr variable LR (default, NULL).
#' @param treatment.names Optional. Names/IDs of the treatment samples, which
#'   must be include in the variable LR (default, NULL).
#' @return The same LR object with a column named "class" added to a GRanges
#'   object from LR (default). Based on the model prediction each DIMP is
#'   labeled as control "CT" or as treatment "TT". If "conf.matrix" is TRUE and
#'   the arguments control.names and treatment.names are provided, then the
#'   overall confusion matrix is returned
#' @examples
#' # Random generation of Hellinger divergence values from a Weibul
#' # distribution model and estimating their tail probabilities.
#' num.points <- 5000
#' set.seed(123)
#' hdiv11 = rweibull(1:num.points, shape = 0.45, scale = 1.2)
#' wprob11 = pweibull(hdiv11, shape = 0.45, scale = 1.2, lower.tail = FALSE)
#' hdiv12 = rweibull(1:num.points, shape = 0.45, scale = 1.2)
#' wprob12 = pweibull(hdiv12, shape = 0.45, scale = 1.2, lower.tail = FALSE)
#' hdiv21 = rweibull(1:num.points, shape = 0.6, scale = 1.02)
#' wprob21 = pweibull(hdiv21, shape = 0.6, scale = 1.02, lower.tail = FALSE)
#' hdiv22 = rweibull(1:num.points, shape = 0.61, scale = 1.02)
#' wprob22 = pweibull(hdiv22, shape = 0.61, scale = 1.02, lower.tail = FALSE)
#' #' Potential signal
#' PS <- GRangesList(
#'       sample11 = makeGRangesFromDataFrame(
#'           data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
#'                      strand = '*', hdiv = hdiv11, wprob = wprob11),
#'           keep.extra.columns = TRUE),
#'           sample12 = makeGRangesFromDataFrame(
#'           data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
#'                      strand = '*', hdiv = hdiv12, wprob = wprob12),
#'           keep.extra.columns = TRUE),
#'       sample21 = makeGRangesFromDataFrame(
#'           data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
#'                      strand = '*', hdiv = hdiv21, wprob = wprob21),
#'           keep.extra.columns = TRUE),
#'           sample22 = makeGRangesFromDataFrame(
#'           data.frame(chr = "chr1", start = 1:num.points, end = 1:num.points,
#'                      strand = '*', hdiv = hdiv22, wprob = wprob22),
#'           keep.extra.columns = TRUE))
#' cutpoint = 5.76
#' DIMPs = selectDIMP(PS, div.col = 1, cutpoint = cutpoint)
#'
#' #' A classification model can be fitted as follow:
#' conf.mat <- evaluateDIMPclass(DIMPs,
#'                               column = c(hdiv = TRUE, TV = FALSE,
#'                                          wprob = FALSE, pos = FALSE),
#'                               interaction = "wprob:hdiv",
#'                               control.names = c("sample11", "sample12"),
#'                               treatment.names = c("sample21", "sample22"))
#' # Now predictions of DIMP for control and treament can be obtained
#' pred = predictDIMPclass(LR = DIMPs, model = conf.mat$model,
#'                         conf.matrix = TRUE,
#'                         control.names = c("sample11", "sample12"),
#'                         treatment.names = c("sample21", "sample22"))
#' pred
#' @importFrom biovizBase flatGrl
#' @export
predictDIMPclass <- function(LR, model, conf.matrix = FALSE,
                             control.names = NULL,
                             treatment.names = NULL) {
  if (conf.matrix && (is.null(control.names) || is.null(treatment.names))) {
    stop(paste0("* if conf.mat = TRUE, then the character vectors for ",
                "control.names and treatment.names must be provided"))
  }

  if (conf.matrix && !is.null(control.names) && !is.null(treatment.names)) {
    sn = names(LR)
    idx.ct = match(control.names, sn)
    idx.tt = match(treatment.names, sn)
    CT = GRangesList(LR[ idx.ct ])
    CT = flatGrl(CT)
    TT = GRangesList(LR[ idx.tt ])
    TT = flatGrl(TT)
    classSet = list(CT = CT, TT = TT)
  } else classSet = LR

  position <- function(gr) {
    chrs = split(gr, seqnames(gr))
    gr = lapply(chrs, function(grc) {
      x = start(grc)
      x.min = min(x)
      x.max = max(x)
      delta =  max(c(x.max - x, 1))
      return((x - x.min)/(delta))
    })
    return(unlist(gr))
  }

  classifier <- function(GR) {
    cn = colnames(mcols(GR))
    vn = c("hdiv", "TV", "wprob")
    vn = cn[is.element(cn, vn)]
    m = mcols(GR[, vn])
    m$pos = position(GR)
    if (is.element("wprob", vn)) {
      vn = sub("wprob", "logP", vn)
      colnames(m) <- c(vn, "pos")
      m$logP = log10(m$logP + 2.2e-308)
    }

    if (inherits(model, "glm")) {
      model = structure(model, class = c("LogisticR", "glm"))
    }
    if (inherits(model, "lda") || inherits(model, "qda")) {
      pred = predict(model, newdata = m )$class
    } else pred = predict(model, newdata = m, type = "class")
    GR$class = pred
    return(GR)
  }
  LR = lapply(classSet, classifier)
  if (!conf.matrix) {
    return(LR)
  } else {
    conf.mat = data.frame(TRUE.class = c(rep("CT", length(CT)),
                                         rep("TT", length(TT))),
                          PRED.class = c(LR$CT$class, LR$TT$class))
    conf.mat = table(conf.mat)
    return(list(conf.mat = conf.mat,
                accuracy = sum(diag(conf.mat)/sum(conf.mat))))}
}

