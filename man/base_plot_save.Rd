% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{base_plot_save}
\alias{base_plot_save}
\title{Generate suitable output graphics device by file suffix.}
\usage{
base_plot_save(saveplot, ...)
}
\arguments{
\item{saveplot}{Save plot to given file "a.pdf", "b.png".}

\item{...}{Additional parameters given to plot output (\code{\link{pdf}}, \code{\link{png}},...) like "width", "height", .etc.}
}
\description{
Generate suitable output graphics device by file suffix.
}
\examples{

base_plot_save("a.pdf")
# will simplify run (pdf("a.pdf))

}
