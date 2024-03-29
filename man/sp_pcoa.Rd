% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_pcoa.R
\name{sp_pcoa}
\alias{sp_pcoa}
\title{Generating pcoa plot}
\usage{
sp_pcoa(
  data,
  metadata,
  input_type = "normalized_OTUtable",
  dissimilarity_index = "bray",
  k = 3,
  top_n = 1,
  statistical_value_type = mad,
  binary_dissimilarity_index = F,
  data_transform = "auto",
  group_variable = NULL,
  color_variable = NULL,
  color_variable_order = NULL,
  shape_variable = NULL,
  shape_variable_order = NULL,
  size_variable = NULL,
  size_variable_order = NULL,
  label_variable = NULL,
  label_variable_order = NULL,
  legend.position = "right",
  draw_ellipse = "auto",
  manual_color_vector = NULL,
  title = NULL,
  label_font_size = NULL,
  debug = FALSE,
  type = "t",
  level = 0.95,
  filename = NULL,
  extra_ggplot2_cmd = NULL,
  check_significance = T,
  check_paired_significance = T,
  coord_fixed = T,
  ...
)
}
\arguments{
\item{data}{Data file (or data.frame) With the first row as the header line and the first column as the row-name.
Columns are separated by one \code{tab}. Each row represents one variable (normally genes, OTU).
Each column represents one sample. The numbers represent gene expression abundance or OTU abundance or
other abundances, and should be normalized.}

\item{metadata}{Metadata file (or data.frame) with sample attributes like group information.
The first column is the same as the first row of value given to parameter \code{data}.
These attributes would be used as \code{color}, \code{size}, \code{shape} variables in the plot.
If not supplied, each sample will be treated as one group.}

\item{input_type}{The input data is OTU table (\code{normalized_OTUtable}) or
a distance matrix (\code{distance_matrix}).}

\item{dissimilarity_index}{Dissimilarity index, partial match to "manhattan", "euclidean",
"canberra", "clark", "bray" (default, meaning Bray–Curtis), "kulczynski", "jaccard", "gower", "altGower",
"morisita", "horn", "mountford", "raup", "binomial", "chao", "cao",
"mahalanobis", "chisq" or "chord".
Gower, Bray–Curtis, Jaccard and Kulczynski indices are good in detecting underlying
ecological gradients (Faith et al. 1987).
Morisita, Horn–Morisita, Binomial, Cao and Chao indices should be able to handle different
sample sizes (Wolda 1981, Krebs 1999, Anderson & Millar 2004),
and Mountford (1962) and Raup-Crick indices for presence–absence data should be able
to handle unknown (and variable) sample sizes.}

\item{k}{the maximum dimension of the space which the data are to be
    represented in; must be in \eqn{\{1, 2, \ldots, n-1\}}{{1, 2, \\dots, n-1}}.}

\item{top_n}{An integer larger than 1 will be used to get top x genes (like top 5000).
A float number less than 1 will be used to get top x fraction genes (like top 0.7 of
all genes).}

\item{statistical_value_type}{Specify the way for statistical computation. Default mad, accept mean, var, sum, median.}

\item{binary_dissimilarity_index}{Perform presence/absence standardization
before computing \code{dissimilarity_index}. Default \code{FALSE}, accept \code{TRUE}.}

\item{data_transform}{Methods for transforming data. Default 'auto'. Accept 'None',
For \code{auto}: If the data values are larger than common abundance class scales
(here is \code{9}), the function performs a Wisconsin double standardization (\code{\link[vegan]{wisconsin}}).
If the values look very large, the function also performs \code{\link{sqrt}} transformation.
For \code{None}: No transformation would be performed.
For \code{total}: Compute relative abundance of OTUs/Genes in each sample.
For \code{hellinger}: square root of \code{method = "total"}.
For \code{scale}: row scale.
For \code{sqrt} and \code{log2}, just as the words.}

\item{group_variable}{The variable for grouping points to generate normal data confidence ellipses. Optional.}

\item{label_variable}{The variable for text used to label points. Optional.
Specially supplying \code{Row.names} would label sample with their names.}

\item{legend.position}{Position of legend, accept top, bottom, left, right, none or c(0.8,0.8).}

\item{draw_ellipse}{Default 'auto' means to enclose points in a polygon if one group with
less than 4 points. If there are more than 4 points for all groups, confidence ellipses would be draw.
Accept \verb{confidence ellipse} to draw confidence ellipses for all conditions even though they would
not be draw.
Accept \code{no} to remove ellipse or other polygens.}

\item{manual_color_vector}{Manually set colors for each geom.
Default NULL, meaning using ggplot2 default.
Colors like c('red', 'blue', '#6181BD') (number of colors not matter) or
a RColorBrewer color set like  "BrBG"     "PiYG"     "PRGn"     "PuOr"
"RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"  "Spectral" "Accent"
"Dark2"    "Paired"   "Pastel1"  "Pastel2"  "Set1"
"Set2"    "Set3"     "Blues"    "BuGn"     "BuPu"
"GnBu"     "Greens"   "Greys"    "Oranges" "OrRd"     "PuBu"
"PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"
"YlGn"    "YlGnBu"   "YlOrBr"   "YlOrRd"
(check http://www.sthda.com/english/wiki/colors-in-r for more).}

\item{title}{Title of picture.}

\item{type}{The type of ellipse.
The default \code{"t"} assumes a multivariate t-distribution, and
\code{"norm"} assumes a multivariate normal distribution.
\code{"euclid"} draws a circle with the radius equal to \code{level},
representing the euclidean distance from the center.
This ellipse probably won't appear circular unless \code{coord_fixed()} is applied.}

\item{level}{The level at which to draw an ellipse,
or, if \code{type="euclid"}, the radius of the circle to be drawn.}

\item{filename}{Output picture to given file.}

\item{extra_ggplot2_cmd}{Extra ggplot2 commands (currently unsupported)}

\item{check_significance}{Check if the centroids and dispersion of the groups
as defined by measure space are equivalent for all groups.}

\item{check_paired_significance}{Paired-check for each two groups.}

\item{coord_fixed}{When True (the default) ensures that one unit on the x-axis is
the same length as one unit on the y-axis.}

\item{...}{Parameters given to \code{sp_ggplot_layout}}
}
\value{
A ggplot2 object
}
\description{
Generating pcoa plot
}
\examples{

## Not run:
## End(Not run)

}
