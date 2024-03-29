% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_raincloud.R
\name{sp_raincloud}
\alias{sp_raincloud}
\title{raincloud}
\usage{
sp_raincloud(
  data,
  melted = TRUE,
  xvariable = NULL,
  yvariable = NULL,
  metadata = NULL,
  ID_var = c(),
  coordinate_flip = TRUE,
  position_nudge_flat_violin_x = 0.3,
  position_nudge_flat_violin_y = 0,
  position_nudge_flat_violin_alpha = 0.8,
  palette_color = "Set2",
  palette_fill = "Set2",
  position_nudge_box_x = 0.25,
  box_fill = "white",
  legend.position = NULL,
  extra_ggplot2_cmd = NULL,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  additional_theme = NULL,
  debug = F,
  ...
)
}
\arguments{
\item{data}{Data file (with header line, the first row is the colname,
tab seperated. Multiple formats are allowed and described above)}

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

\item{metadata}{Giving a metadata file with format specified in example
to tell the group information for each sample.}

\item{ID_var}{Other columns one want to treat as ID variable columns
except the one given to \code{xvariable}.}

\item{coordinate_flip}{Rotate the plot from vertical to horizontal.
Usefull for plots with many values or very long labels at X-axis}

\item{position_nudge_flat_violin_x}{The violin moves on the X-axis. Default 0.3.}

\item{position_nudge_flat_violin_y}{The violin moves on the Y-axis. Default 0.}

\item{position_nudge_flat_violin_alpha}{The violin transparency.}

\item{palette_color}{The violin palette.}

\item{palette_fill}{The point palette.}

\item{position_nudge_box_x}{The box moves on the X-axis. Default 0.25.}

\item{box_fill}{Box color.}

\item{...}{Parametes given to \code{sp_ggplot_layout}}
}
\value{
A ggplot2 object
}
\description{
raincloud
}
\examples{
NULL

}
