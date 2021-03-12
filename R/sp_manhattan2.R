


#' Generating manhattan plot
#'
#' @param data Data frame or data file (with header line, the first column will not be treated as the rowname, tab seperated).
#' @param ID_var Name of ID column.
#' @param FDR_var Name of FDR column.
#' @param grp_var Group variable for ordering points.
#' @param sig_col Significant column.Optional, in sample data significant.
#' @param status_col_level Changing the order of status column values. Normally, the unique values of status column would be sorted alphabetically.
#' @param pvalue Set the filter threshold for defining significance. Default "0.05".
#' @param shape_col Column for points shapes.
#' @param color_var Column for color points. Default the same as group variable.
#' @param alpha Transparent alpha value.Default 0.4, Accept a float from 0(transparent), 1(opaque).
#' @param grp_var_order  Group variable order, like "'Rhizobiales', 'Actinomycetales'".
#' @param color_var_order Color variable order. Default same as `grp_var`.
#' @param shape_col_order Levels for shapes, like "'Sig','nonSig'".
#' @param log10_transform_fdr Get `-log10(FDR)` for column given to `FDR_var`. Default FALSE, accept TRUE.
#' @param point_label_var Name of columns containing labels for points.
#' @param x_label Xlab label.Default NULL.
#' @param y_label Ylab label.Default "Negative log10 FDR".
#' @inheritParams sp_ggplot_add_vline_hline
#' @inheritParams sp_ggplot_layout
#' @param ...
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#'
#' ## Not run:
#' manhattan_data = "manhattan.data"
#'
#' sp_manhattan2_plot(data=manhattan_data, ID_var='ID', FDR_var='FDR', title="test1", point_size=2, point_label_var = "Labels")
#' ## End(Not run)

sp_manhattan2_plot <- function(data,
                               ID_var = "ID",
                               FDR_var = "FDR",
                               grp_var = 'Phylum',
                               sig_col = "level",
                               status_col_level = c(),
                               pvalue = 0.05,
                               shape_col = NULL,
                               color_var = NULL,
                               alpha = 0.4,
                               grp_var_order = c(),
                               color_var_order = c(),
                               shape_col_order = c(),
                               log10_transform_fdr = TRUE,
                               filename = NULL,
                               point_label_var = 'CTctctCT',
                               title='',
                               x_label = NULL,
                               y_label = "Negative log10 FDR",
                               legend.position = "right",
                               point_size = 'ct___',
                               ...) {
  options(warn = -1)
  options(scipen = 999)

  if (class(data) == "character") {
    data <- sp_readTable(data, row.names = NULL)
  }

  if (is.null(sig_col) && pvalue != '') {
    sig_col = 'Significant'
    self_compute_status = TRUE
  } else {
    self_compute_status = FALSE
  }

  if (is.null(shape_col)) {
    shape_col = sig_col
  }

  if (is.null(color_var)) {
    color_var = grp_var
  }

  if (!self_compute_status) {
    if (length(status_col_level)) {
      sig_level <- status_col_level
    } else {
      sig_level = c()
    }
  } else {
    sig_level <- c("Sig", "Unsig")
    data[[sig_col]] <-
      ifelse(data[[FDR_var]] <= pvalue, "Sig", "Unsig")
    data$alpha <-
      ifelse(data[[FDR_var]] <= pvalue, alpha, alpha / 2)
  }

  if (length(sig_level) > 1) {
    data[[sig_col]] <-
      factor(data[[sig_col]], levels = sig_level, ordered = T)
  }

  if (is.numeric(data[[shape_col]])) {
    stop("Shape variable must not be numerical columns.")
  }


  if (length(grp_var_order) > 1) {
    data[[grp_var]] <-
      factor(data[[grp_var]], levels = grp_var_order, ordered = T)
  }

  if (length(color_var_order) > 1) {
    data[[color_var]] <-
      factor(data[[color_var]], levels = color_var_order, ordered = T)
  }

  if (length(shape_col_order) > 1) {
    data[[shape_col]] <-
      factor(data[[shape_col]], levels = shape_col_order, ordered = T)
  }

  if (is.numeric(data[[FDR_var]]) && log10_transform_fdr) {
    data[[FDR_var]] <- (-1) * log10(data[[FDR_var]])
  } else {
    stop(
      "Y axis variable must be numeric , normally it should be p-value or adjusted p-value or similar columns."
    )
  }

  data <- data[order(data[[grp_var]], data[[ID_var]]), ]

  data[[ID_var]] <- paste(data[[grp_var]], data[[ID_var]], sep = "")

  data2 <- data.frame(grp = data[[grp_var]], num = 1:nrow(data))
  data2_mean <- data2 %>% group_by(grp) %>% summarize_all(median)
  data2_mean <- data2_mean[order(data2_mean$num), ]

  data[[ID_var]] <- data2[, 2]


  if (color_var == grp_var) {
    show_color_legend = 'none'
  } else {
    show_color_legend = 'legend'
  }

  if (shape_col == grp_var) {
    show_shape_legend = 'none'
  } else {
    show_shape_legend = 'legend'
  }

  ID_var_en = sym(ID_var)
  FDR_var_en = sym(FDR_var)
  color_var_en = sym(color_var)
  shape_col_en = sym(shape_col)
  p <-
    ggplot(data = data, aes(x = !!ID_var_en, y = !!FDR_var_en)) + theme_classic()

  p <- p + geom_point(alpha = alpha) + labs(title = title) +
    aes(colour = !!color_var_en,
        shape = !!shape_col_en) +
    guides(colour = show_color_legend, shape = show_shape_legend)

  if (point_size != "ct___") {
    p <- p + aes(size = point_size)
    if (is.numeric(data[[point_size]])) {
      p <- p + scale_size_continuous(range = c(0.1, 2))
    }
  }


  if (point_label_var != "CTctctCT") {
    data.l <-
      data[data[[point_label_var]] != "-" &
             data[[point_label_var]] != "", ]
    p <-
      p + geom_text(
        data = data.l,
        aes(
          x = data.l[[ID_var]],
          y = data.l[[FDR_var]],
          label = data.l[[point_label_var]]
        ),
        colour = "black",
        show.legend = F
      )
  }

  p <-
    p + scale_x_continuous(breaks = data2_mean$num, labels = data2_mean$grp) +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ))


  p <- sp_ggplot_add_vline_hline(p,
                                 yintercept = (-1) * log10(pvalue))


  p <- sp_ggplot_layout(
    p,
    xtics_angle = 0,
    coordinate_flip= FALSE,
    legend.position = legend.position,
    extra_ggplot2_cmd = ,
    filename = filename,
    title = title,
    x_label = x_label,
    y_label = y_label,
    ...)

  p

}
