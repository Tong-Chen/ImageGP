#!/usr/bin/env Rscript




data = "system_file/histogram/histogram_demo2.txt"
melted = FALSE
xvariable = "NULL"
color_variable = "Cell"
color_variable_order = "NULL"
group_variable = "NULL"
group_variable_order = "NULL"
value_scale = 1
facet_variable = "NULL"
facet_variable_order = "NULL"
maximum_allowed = Inf
yaxis_statistics = "count/sum(count)"
alpha = 0.4
plot_type = "line"
binwidth = NULL
hist_bar_position = "identity"
yaxis_scale_mode = "NULL"
y_add = 0
fill_area = FALSE
line_size = 1
add_mean_value_vline = FALSE
facet_scales = "fixed"
facet_ncol = NULL
facet_nrow = NULL
xtics = TRUE
ytics = TRUE
legend.position = "right"
xtics_angle = 0
manual_xtics_value = "NULL"
manual_xtics_pos = "NULL"
custom_vline_x_position = "NULL"
custom_vline_anno = "NULL"
manual_color_vector = "Accent"
x_label = "NULL"
y_label = "NULL"
title = "NULL"
width = 15
height = 15
outputprefix = "system_file/histogram/histogram_demo2.txt"
outputpictype = "pdf"
saveppt = FALSE
savehtml = FALSE


# For all users
# devtools::install_github("Tong-Chen/ImageGP")
# For chinese users
# devtools::install_git("https://gitee.com/ct5869/ImageGP.git")
library(ImageGP)
library(ggplot2)
library(reshape2)
library(grid)
library(dplyr)
library(htmlwidgets)
library(plotly)

if (data == "") {
  script = sub(".*=", "", commandArgs()[4])
  #print(script)
  system(paste(script, "-h"))
  stop("At least -f is required!")
}



if (outputprefix == "") {
  outputprefix = data
}

filename = paste0(outputprefix, '.histogram.', outputpictype)


color_variable_order = sp_string2vector(color_variable_order)
group_variable_order = sp_string2vector(group_variable_order)
facet_variable_order = sp_string2vector(facet_variable_order)
manual_color_vector = sp_string2vector(manual_color_vector)
custom_vline_x_position = sp_string2vector(custom_vline_x_position)
manual_xtics_pos = sp_string2vector(manual_xtics_pos)
manual_xtics_value = sp_string2vector(manual_xtics_value)
custom_vline_anno = sp_string2vector(custom_vline_anno)

cat(sp_current_time(), "Starting...\n")

sp_histogram(
  data = data,
  melted = melted,
  xvariable = xvariable,
  color_variable = color_variable,
  color_variable_order = color_variable_order,
  group_variable = group_variable,
  group_variable_order = group_variable_order,
  value_scale = value_scale,
  facet_variable = facet_variable,
  facet_variable_order = facet_variable_order,
  maximum_allowed = maximum_allowed,
  yaxis_statistics = yaxis_statistics,
  plot_type = plot_type,
  binwidth = binwidth,
  hist_bar_position = hist_bar_position,
  alpha = alpha,
  yaxis_scale_mode = yaxis_scale_mode,
  y_add = y_add,
  fill_area = fill_area,
  line_size = line_size,
  add_mean_value_vline = add_mean_value_vline,
  manual_color_vector = manual_color_vector,
  facet_scales = facet_scales,
  facet_ncol = facet_ncol,
  facet_nrow = facet_nrow,
  xtics = xtics,
  ytics = ytics,
  width = width,
  height = height,
  legend.position = legend.position,
  xtics_angle = xtics_angle,
  manual_xtics_value = manual_xtics_value,
  manual_xtics_pos = manual_xtics_pos,
  custom_vline_x_position = custom_vline_x_position,
  custom_vline_anno = custom_vline_anno,
  x_label = x_label,
  y_label = y_label,
  title = title,
  filename = filename,
  saveppt = saveppt,
  savehtml = savehtml
)
cat(sp_current_time(), "Success.\n")


