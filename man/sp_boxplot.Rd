% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_boxplot.R
\name{sp_boxplot}
\alias{sp_boxplot}
\title{Generating box plot}
\usage{
sp_boxplot(
  data,
  melted = FALSE,
  xvariable = NULL,
  yvariable = NULL,
  legend_variable = NULL,
  statistics = FALSE,
  xtics_angle = 0,
  legend_variable_order = NULL,
  legend_variable_cut = NULL,
  xvariable_order = NULL,
  xvariable_cut = NULL,
  xvariable_cut_order = NULL,
  group_variable_for_line = NULL,
  group_variable_order_for_line = NULL,
  y_add = 0,
  yaxis_scale_mode = NULL,
  notch = FALSE,
  par = NULL,
  outlier = FALSE,
  statistical_method = "aov",
  statistical_threshold_for_letters = 0.05,
  out_scale = 1.05,
  legend.position = "right",
  manual_color_vector = NULL,
  violin = FALSE,
  violin_nb = FALSE,
  scale_violin = "width",
  ID_var = c(),
  jitter = FALSE,
  jitter_bp = FALSE,
  dotplot = FALSE,
  dotplot_bp = FALSE,
  colormodel = "srgb",
  coordinate_flip = FALSE,
  facet_variable = NULL,
  facet_variable_order = NULL,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  facet_nrow = NULL,
  facet_ncol = NULL,
  facet_singlecell_style = F,
  facet_scales = "fixed",
  metadata = NULL,
  debug = F,
  filename = NULL,
  extra_ggplot2_cmd = NULL,
  dotsize = 2,
  ...
)
}
\arguments{
\item{data}{Data file (with header line, the first row is the colname,
tab separated. Multiple formats are allowed and described above)}

\item{melted}{When TRUE, meaning a long format matrix is supplied to \code{data}.
function will skip preprocess. Default FALSE.}

\item{xvariable}{The column represents the x-axis values. For unmelted data, the program
will use first column as x-variable. If one want to use first row of unmelted data
as x-variable, please specify \code{variable} here (which is an inner name).
Or if one want to use other columns in \code{metadata}.}

\item{yvariable}{The column represents the digital values.
For unmelted data, the program
will use \code{value} as y-variable (which is an inner name).
This parameter can only be set when \code{melted} is TRUE.}

\item{legend_variable}{The column represents the legend information.
Default \code{xvariable} if not specified.}

\item{xtics_angle}{Rotation angle for a-axis. Default 0.}

\item{legend_variable_order}{Levels for legend variable.
Default data order, accept a vector like c('TP16','TP22','TP23') for \code{legend_variable} column.}

\item{legend_variable_cut}{Self-define intervals for legend variable when
values in \code{legend_variable} column is continuous numbers.}

\item{xvariable_order}{xvariable_order Levels for x-axis variable. Default data order, accept input like c('g','a','j','x','s','c','o','u') for Set column.}

\item{xvariable_cut}{xvariable_cut Self-define intervals for x-axis variable.}

\item{xvariable_cut_order}{Specify a list of names for self-define intervals for x-axis variable.}

\item{group_variable_for_line}{Specify the group of points to line together (one column name).}

\item{group_variable_order_for_line}{Levels for group variable for lines.}

\item{y_add}{A number to add if log scale is used.
Default 0 meaning the minimum non-zero value would be used.}

\item{yaxis_scale_mode}{Give the following \code{scale_y_log10()},
\code{coord_trans(y="log10")}, or other legal command for ggplot2 or
simply \code{log2} to set the scale way.}

\item{notch}{Using notch (sand clock shape) or not. Default FALSE.}

\item{outlier}{Exclude outliers. Exclude outliers or not, default `FALSE`` means keeping outliers.}

\item{statistical_method}{Statistical method. For two groups, default <t.test>. For more than two groups, default \if{html}{\out{<aov>}}.}

\item{statistical_threshold_for_letters}{Threshold for treating as significance, default 0.05.}

\item{out_scale}{The scales for one want to set to exclude outliers.
Default 1.05. No recommend to change unless you know what you are doing.}

\item{legend.position}{Position of legend, accept top, bottom, left, right, none or c(0.8,0.8).}

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

\item{violin}{Do violin plot plus inner boxplot.}

\item{violin_nb}{Do violin plot without inner boxplot.}

\item{scale_violin}{The value given to scale for violin plot.
if "area", all violins have the same area (before trimming the tails).
If "count", areas are scaled proportionally to the number of observations.
If "width", all violins have the same maximum width. 'equal' is also accepted.}

\item{ID_var}{Other columns one want to treat as ID variable columns
except the one given to \code{xvariable}.}

\item{jitter}{Do jitter plot instead of boxplot.}

\item{jitter_bp}{Do jitter plot overlay with violin plot or boxplot or both.}

\item{dotplot}{Do dotplot plot instead of boxplot.}

\item{dotplot_bp}{Do dotplot plot overlay with violin plot or boxplot or both.}

\item{coordinate_flip}{Rotate the plot from vertical to horizontal.
Usefull for plots with many values or very long labels at X-axis}

\item{facet_variable}{Wrap plots by given column (one of column names should be specified).
This is used to put multiple plot in one picture.}

\item{facet_variable_order}{The levels of wrapping to set the order of each group.}

\item{x_label}{Xlab label.}

\item{y_label}{Ylab label.}

\item{title}{Title of picture.}

\item{facet_nrow}{The number of rows one want when \code{facet} is used. Default NULL.}

\item{facet_ncol}{The number of columns one want when \code{facet} is used. Default NULL.}

\item{facet_singlecell_style}{Use specified style for showing single cell gene expression profile. Default FALSE.}

\item{facet_scales}{Paramter for scales for facet. Default \code{fixed} meaning each inner graph
use same scale (x,y range), \code{free} (variable x, y ranges for each sub-plot),
\code{free_x} (variable x ranges for each sub-plot), \code{free_y} (variable y ranges for each sub-plot).}

\item{metadata}{Giving a metadata file with format specified in example
to tell the group information for each sample.}

\item{filename}{Output picture to given file.}

\item{extra_ggplot2_cmd}{Extra ggplot2 commands (currently unsupported)}

\item{...}{Parametes given to \code{sp_ggplot_layout}}

\item{The}{size of dots. \link{Default 2}}
}
\value{
A ggplot2 object
}
\description{
\code{metadata}
}
\examples{

box_test_data <- data.frame(ID = letters[1:4],
Gene = letters[c(8,8,9,9,10,10,11,11)], Expr = runif(16))
sp_boxplot(data = box_test_data, xvariable = "ID", value = "Expr", variable ="Gene")

## Not run:
box_data = "box.data"

sp_boxplot(data = box_data, xvariable = "Gene", value = "Expr", variable="Group")
## End(Not run)

}
