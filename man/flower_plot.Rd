% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowerplot.R
\name{flower_plot}
\alias{flower_plot}
\title{Flower plot could be treated as one kind of venn diagram but only showing common items
like OTUs or genes among all groups and total items (or total items excluding common
items) for each group.}
\usage{
flower_plot(
  input,
  sep = "\\t",
  row.names = NULL,
  header = T,
  quote = "",
  comment = "",
  check.names = F,
  item_variable = NULL,
  set_variable = NULL,
  start = 90,
  a = 0.5,
  b = 2,
  r = 1,
  group_color = "Spectral",
  group_color_alpha = 0.6,
  label_total_num_items = TRUE,
  saveplot = NULL,
  label = "core",
  common_color = "white",
  common_color_alpha = 1,
  saveppt = FALSE,
  ...
)
}
\arguments{
\item{input}{Input data file (first line as header line, the first column is the name of
genes or OTUs or otehr things one wants to compare, the second column is the group name which genes belong to,
tab seperated)

\if{html}{\out{<div class="sourceCode">}}\preformatted{Gene    Sample
g1    Set1
a1    Set3
b4    Set1
.
.
c1    Set3
}\if{html}{\out{</div>}}}

\item{start}{Start position of first ellipse. Default 90 represents starting from 0 clock.}

\item{a, b}{
Vectors, radii of the ellypses along the two axes in user units.
}

\item{r}{Set the size of the center circle.}

\item{group_color}{Set the color of the petal ellipse (each group), with input format，like：c('#6181BD4E','#F348004E','#64A10E4E'...) or
supply a RColorBrewer color set like "Set1", "Set2", "Set3", "YlOrRd"
(check http://www.sthda.com/english/wiki/colors-in-r for more).}

\item{group_color_alpha}{The transparency of each ellipse color. Default 0.6.}

\item{label_total_num_items}{Label total number of items in for each group (when True) or label number of items in each group
after substracting numbe of common items.}

\item{saveplot}{Save plot to given file "a.pdf", "b.png".}

\item{label}{The name of the center circle.}

\item{common_color}{The color of the center circle. Default "white".}

\item{...}{Parameters givento \code{\link{base_plot_save}}}

\item{common_col_alpha}{The transparency of common circle color. Default 0.6.}
}
\value{
An image
}
\description{
Flower plot could be treated as one kind of venn diagram but only showing common items
like OTUs or genes among all groups and total items (or total items excluding common
items) for each group.
}
\examples{

flowerinput <- "test.file"
flower(flowerinput)

}
