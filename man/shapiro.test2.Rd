% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{shapiro.test2}
\alias{shapiro.test2}
\title{An inner version of \code{shapiro.test} with two modifications.}
\usage{
shapiro.test2(data, threshold = 1000)
}
\arguments{
\item{data}{A vector of observation values.}

\item{threshold}{A number to define \strong{large} sample size. Default \code{1000}. For data with
more than given \code{threshold} samples, always assume they past normality test,}
}
\value{
A logical value
}
\description{
\itemize{
\item Always return FALSE for sample size less than 3.
\item Always return TRUE for sample size larger than given threshold (default 1000)
}
}
\examples{

shapiro.test2(c(1,2,3,4))

}
