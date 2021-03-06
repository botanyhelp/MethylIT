% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmpClusters.R
\name{dmpClusters}
\alias{dmpClusters}
\title{DMP clustering}
\usage{
dmpClusters(
  GR,
  maxDist = 3,
  minNumDMPs = 1,
  chromosomes = NULL,
  num.cores = 1L,
  tasks = 0L,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{GR}{An object from \strong{\emph{'pDMP'}} class, which is returned by
\code{\link{selectDIMP}} function.}

\item{maxDist}{maximum distance at which two reported bases sites from the
same cluster can be separated. Default: \eqn{maxDist = 3}.}

\item{minNumDMPs}{minimum number of marked bases inside of each cluster.
Default: \eqn{minNumDMPs = 1}.}

\item{chromosomes}{vector of characters labeling the chromosomes included in
the analysis. Default: chromosomes = NULL (all chromosomes are included).}

\item{num.cores, tasks}{integer(1). The number of cores to use, i.e. at most
how many child processes will be run simultaneously (see
\code{\link[BiocParallel]{bplapply}} function from BiocParallel
package).The number of tasks per job. value must be a scalar integer >=
0L (see MulticoreParam from BiocParallel package).}

\item{verbose}{if TRUE, prints the function log to stdout.}

\item{...}{Further parameters for uniqueGRanges function.}
}
\value{
A GRanges object with the numbers of positions inside each cluster,
where DMPs were reported in at least one of the samples.
}
\description{
Given a 'pDMP' object carrying DMPs obtained in Methyl-IT
downstream analysis, function \strong{\emph{"dmpClusters"}} build
clusters of DMPs, which can be further tested to identify differentially
methylated regions (DMRs) with \code{\link{countTest2}} function.
}
\details{
DNA base coordinates provided in the \strong{\emph{GR}} object are
used to build a binary string of zeros and ones . The binary string of zeros
and ones is used in the detection of cluster of DMPs. Postulating that a
reported DNA base is found at the beginning and the end of a cluster,
genome-wide screening can be performed where two consecutive targeting base
positions are separated by less than a given threshold
\strong{\emph{maxDist}} bases (1). DMPs from all the samples included in
\strong{\emph{GR}} are considered to build the clusters.

The number of DMPs reported in each cluster corresponds to the numbers of
positions inside the cluster where DMPs were reported in at least one
of the samples.
}
\examples{
## Get a dataset of potential signals and the estimated cutpoint from the
## package
data(PS, cutpoint)

## The estimated cutpoints are used to discriminate signals from the noise.
## That is, DMPs are selected using the cupoints
dmps <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)

## Build clusters of DMPs
x1 = dmpClusters(GR = dmps, maxDist = 7, minNumDMPs = 6, num.cores=2L,
                verbose = FALSE)
x1

## ==== Setting up the experiment design to test for DMRs ===
nams <- names(dmps)
dmps_at_clusters <- getDIMPatGenes(GR = dmps, GENES = x1,
                                    ignore.strand = TRUE)
dmps_at_clusters <- uniqueGRanges(dmps_at_clusters, columns = 2L,
                                    ignore.strand = TRUE, type = "equal",
                                    verbose = FALSE)
colnames(mcols(dmps_at_clusters)) <- nams

colData <- data.frame(condition = factor(c("CT", "CT", "CT",
                                           "TT", "TT", "TT"),
                                         levels = c("CT", "TT")),
                    nams,
                    row.names = 2)

## Build a RangedGlmDataSet object
ds <- glmDataSet(GR = dmps_at_clusters, colData = colData)

## ================ Testing for DMRs ===========
dmrs <- countTest2(DS = ds, num.cores = 4L,
                   minCountPerIndv = 4,
                   maxGrpCV = c(1, 1),
                   Minlog2FC = 0.5,
                   CountPerBp = 0.001,
                   test = "LRT",
                   verbose = TRUE)
dmrs
}
\references{
\enumerate{
\item Sanchez R, Mackenzie SA (2016) Information Thermodynamics of
Cytosine DNA Methylation. PLOS ONE 11(3): e0150427.
https://doi.org/10.1371/journal.pone.0150427
}
}
\author{
Robersy Sanchez (\url{https://genomaths.com}).
}
