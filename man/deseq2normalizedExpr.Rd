% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcriptome.R
\name{deseq2normalizedExpr}
\alias{deseq2normalizedExpr}
\title{To output normalized results to files named by given "output_prefix",
and also return a list containing normalized counts for downstream analysis.}
\usage{
deseq2normalizedExpr(
  dds,
  output_prefix = "ehbio",
  rlog = T,
  vst = F,
  savemat = T,
  design = NULL
)
}
\arguments{
\item{dds}{\code{\link{salmon2deseq}} or \code{\link{readscount2deseq}} or
\code{\link[DESeq2]{DESeq}} generated dataset.}

\item{output_prefix}{A string, will be used as output file name prefix.}

\item{rlog}{Get "rlog" transformed value for downstream correlation like analysis.}

\item{vst}{Get "vst" transformed value for downstream correlation like analysis. Normally faster than "rlog".}

\item{savemat}{Save normalized and rlog/vst matrix to file. Default T.
The file would be named like \code{output_prefix.DESeq2.normalized.xls},
\code{output_prefix.DESeq2.normalized.rlog.xls}.}
}
\value{
A list containing normalized expression values and/or rlog, vst transformed normalized expression values.
}
\description{
To output normalized results to files named by given "output_prefix",
and also return a list containing normalized counts for downstream analysis.
}
\examples{

nomrexpr <- deseq2normalizedExpr(dds)

}
