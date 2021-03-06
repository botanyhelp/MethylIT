% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poolFromGRlist.R
\name{poolFromGRlist}
\alias{poolFromGRlist}
\title{Methylation pool from a list of GRanges objects with methylation read
counts}
\usage{
poolFromGRlist(
  LR,
  stat = c("mean", "median", "jackmean", "sum"),
  num.cores = 1,
  tasks = 0L,
  prob = FALSE,
  column = 1L,
  jstat = c("sum", "mean", "median"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{LR}{list of GRanges objects to build a virtual individual (methylation
pool)}

\item{stat}{statistic used to estimate the methylation pool: row "mean", row
"median", row "sum", or Jacknife row mean ("jackmean") of methylated and
unmethylated read counts across individuals. Notice that, for only two
samples, "jackmean" makes not sense. Since the centrality statistics are
sensitive to extreme values, stat = 'sum' is an attractive option.
However, in this last case, a further correction for the minimum coverage
for the reference sample must be taken into account in a furhter
estimation of the Hellinger divergence of methylation levels, which is
explained in the detail section from the help of function
\code{\link{estimateDivergence}}. A conservative option is "mean", which
will return the group centroid.}

\item{num.cores}{The number of cores to use, i.e. at most how many child
processes will be run simultaneously (see bplapply function from
BiocParallel package).}

\item{tasks}{integer(1). The number of tasks per job. Value must be a scalar
integer >= 0L. In this documentation a job is defined as a single call
to a function, such as bplapply, bpmapply etc. A task is the division of
the X argument into chunks. When tasks == 0 (default), X is divided as
evenly as possible over the number of workers (see MulticoreParam from
BiocParallel package).}

\item{prob}{Logic. Whether the variable for pooling is between 0 and 1 (a
probability), e.g., methylation levels. If TRUE, then Fisher's
transformation is applied, the row mean is computed for each cytosine
site and returned in the original measurement scale between 0 and 1 by
using the inverse of Fisher's transformation.}

\item{column}{If prob == TRUE, then the 'column' from the LR metacolumns
where the prob values are found must be provided. Otherwise, column = 1L.}

\item{jstat}{If stat = "jackmean", then any of the 'stat' possible values:
"sum", "mean", or "median" can be used to compute, for each cytosine
site, the Jacknife vector of the selected statistics and then to
compute the corresponding mean. Default is jstat = "sum".}

\item{verbose}{If TRUE, prints the function log to stdout}

\item{...}{Additional parameters for 'uniqueGRanges' function.}
}
\value{
A GRanges object
}
\description{
This function will build a GRanges methylation pool from a list
of GRanges objects
}
\details{
The list of GRanges objects (LR) provided to build a virtual
methylome should be an output of the function 'readCounts2GRangesList' or
at least each GRanges must have the columns named "mC" and "uC", for the
read counts of methylated and unmethylated cytosines, respectively.
}
\examples{
gr1 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = '*', mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(
    data.frame(chr = "chr1", start = 11:15, end = 11:15,
               strand = '*', mC = 1, uC = 1:5),
    keep.extra.columns = TRUE)

answer <- poolFromGRlist(list(gr1, gr2), stat = 'mean', verbose = FALSE)

}
