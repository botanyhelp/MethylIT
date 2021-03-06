% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateECDF.R
\name{estimateECDF}
\alias{estimateECDF}
\title{A variant of Empirical Cumulative Distribution Function "ecdf"}
\usage{
estimateECDF(x, npoints = 10)
}
\arguments{
\item{x}{numeric vector}

\item{npoints}{minimum number of non-missing values}
}
\value{
ecdf of numeric vector
}
\description{
This function is used to reduce the number of points used to
build a ecdf of an experimental variable. When a variable has a very
high amount of experimental values (several millions) could be
computationally time consuming to perform a good-of-fit test and a
non-linear regression estimation for a theoretical CDF based in such a
big number of values.
}
\details{
The histogram cell midpoints values are used to build estimateECDF.
}
\examples{
    x <- sample(1:500, 50, replace=TRUE)
    estimateECDF(x, npoints = 15)

}
\keyword{internal}
