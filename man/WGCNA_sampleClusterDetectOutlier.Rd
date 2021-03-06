% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WGCNA.R
\name{WGCNA_sampleClusterDetectOutlier}
\alias{WGCNA_sampleClusterDetectOutlier}
\title{Sample cluster and outlier detection}
\usage{
WGCNA_sampleClusterDetectOutlier(
  wgcnaL,
  thresholdZ.k = -2.5,
  saveplot = NULL,
  removeOutlier = F,
  traitColors = NULL,
  ...
)
}
\arguments{
\item{wgcnaL}{A matrix or an object return by \code{WGCNA_readindata}. A transformed gene expression matrix normally output by \code{WGCNA_dataFilter}.
Samples x Genes.}

\item{thresholdZ.k}{Threshold for defining outliers. First compute the overall
corelation of one sample to other samples. Then do Z-score transfer for all
correlation values. The samples with corelation values less than given value
would be treated as outliers.
Default -2.5 meaning -2.5 std.}

\item{saveplot}{Save plot to given file "a.pdf", "b.png".}

\item{removeOutlier}{Remove outlier samples. Normally this should be only performed if
no suitable soft power can be found.}

\item{traitColors}{Sample attributes data frame transferred by
\code{\link[WGCNA]{numbers2colors}} or generated in \code{\link{WGCNA_readindata}}.}

\item{...}{Additional parameters given to plot output (\code{\link{pdf}}, \code{\link{png}},...) like "width", "height", .etc.}
}
\value{
A data frame.
}
\description{
Sample cluster and outlier detection
}
\examples{

df = generateAbundanceDF(nSample=30, nGrp=3)
datExpr <- WGCNA_dataFilter(df)
datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)


exprMat <- "test.file"
wgcnaL <- WGCNA_readindata(exprMat)

traitData <- 'trait.file'
wgcnaL <- WGCNA_readindata(exprMat, traitData)
datExpr <- wgcnaL$datExpr
WGCNA_dataCheck(datExpr)
datExpr <- WGCNA_dataFilter(datExpr)
datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)

}
