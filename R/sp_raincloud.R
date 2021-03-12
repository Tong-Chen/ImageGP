# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r

# GeomFlatViolin函数的定义见https://github.com/EasyChart/Beautiful-Visualization-with-R
# https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R

"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}

geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }

#' @inheritParams ggplot2::ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)

      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        group_by(group) %>%
        mutate(
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )

    },

    draw_group = function(data, panel_scales, coord) {
      # Find the points for the line to go all the way around
      data <- transform(data,
                        xminv = x,
                        xmaxv = x + violinwidth * (xmax - x))

      # Make sure it's sorted properly to draw the outline
      newdata <-
        rbind(plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y))

      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])

      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },

    draw_key = draw_key_polygon,

    default_aes = aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),

    required_aes = c("x", "y")
  )


### Example:
# ggplot(diamonds, aes(cut, carat)) +
#     geom_flat_violin() +
#     coord_flip()



#' raincloud
#'
#' @param data Data file (with header line, the first row is the colname,
#' tab seperated. Multiple formats are allowed and described above)
#' @param melted When TRUE, meaning a long format matrix is supplied to `data`.
#' function will skip preprocess. Default FALSE.
#' @param metadata Giving a metadata file with format specified in example
#' to tell the group information for each sample.
#' @param xvariable The column represents the x-axis values. For unmelted data, the program
#' will use first column as x-variable. If one want to use first row of unmelted data
#' as x-variable, please specify `variable` here (which is an inner name).
#' Or if one want to use other columns in `metadata`.
#' @param yvariable The column represents the digital values.
#' For unmelted data, the program
#' will use `value` as y-variable (which is an inner name).
#' This parameter can only be set when `melted` is TRUE.
#' @param ID_var Other columns one want to treat as ID variable columns
#' except the one given to `xvariable`.
#' @param coordinate_flip Rotate the plot from vertical to horizontal.
#' Usefull for plots with many values or very long labels at X-axis
#' @param position_nudge_flat_violin_x The violin moves on the X-axis. Default 0.3.
#' @param position_nudge_flat_violin_y The violin moves on the Y-axis. Default 0.
#' @param position_nudge_flat_violin_alpha The violin transparency.
#' @param palette_color The violin palette.
#' @param palette_fill The point palette.
#' @param position_nudge_box_x  The box moves on the X-axis. Default 0.25.
#' @param box_fill Box color.
#' @param ... Parametes given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
sp_raincloud <- function (data,
                          melted = TRUE ,
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
                          ...) {
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if (!melted) {
    if (sp.is.null(yvariable)) {
      yvariable = "value"
    }
  }

  if (class(data) == "character") {
    data <- sp_readTable(data, row.names = NULL)
    if (!melted) {
      first_column_variable <- colnames(data)[1]
      if (sp.is.null(xvariable)) {
        xvariable = first_column_variable
      }
      data <-
        melt(data, id.vars = c(ID_var, first_column_variable))
    }
  } else if (class(data) == "data.frame") {
    if (!melted) {
      if (sp.is.null(xvariable)) {
        xvariable = colnames(data)[1]
      }
    }
  } else if (class(data) != "data.frame") {
    stop("Unknown input format for `data` parameter.")
  }

  if (sp.is.null(xvariable) || sp.is.null(yvariable)) {
    stop('xvariable or yvariable must be specified!')
  }

  # print(data)
  if (!sp.is.null(metadata)) {
    if (class(metadata) == "character") {
      metadata <- sp_readTable(metadata, row.names = NULL)
    } else if (class(metadata) != "data.frame") {
      stop("Unknown input format for `metadata` parameter.")
    }
    # return(list(data=data, metadata=metadata))
    matched_column <-
      get_matched_columns_based_on_value(data, metadata,
                                         only_allow_one_match =
                                           T)

    # return(list(data=data, metadata=metadata, matched_column=matched_column))
    data <-
      merge(data, metadata, by.x = matched_column[1], by.y = matched_column[2])
  }

  data_colnames <- colnames(data)

  if (!(xvariable %in% data_colnames &&
        yvariable %in% data_colnames)) {
    stop(paste(xvariable, 'or', yvariable, 'must be column names of data!'))
  }


  xvariable_en = sym(xvariable)
  yvariable_en = sym(yvariable)

  theme_boxplot <- theme(
    text = element_text(size = 10),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(vjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    plot.title = element_text(
      lineheight = .8,
      face = "bold",
      size = 16
    ),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(
      colour = 'black',
      size = 0.5,
      linetype = 'solid'
    ),
    axis.line.y = element_line(
      colour = 'black',
      size = 0.5,
      linetype = 'solid'
    )
  )

  p <-
    ggplot(data,
           aes(!!xvariable_en, !!yvariable_en, fill = !!xvariable_en)) +
    geom_flat_violin(
      position = position_nudge(x = position_nudge_flat_violin_x , y = position_nudge_flat_violin_y),
      alpha = position_nudge_flat_violin_alpha
    ) +
    scale_color_brewer(palette = palette_color) +
    scale_fill_brewer(palette = palette_fill) +
    geom_jitter(aes(color = !!xvariable_en), width = 0.1) +
    geom_boxplot(
      width = .1,
      position = position_nudge(x = position_nudge_box_x),
      fill = box_fill,
      size = 0.5
    ) +
    theme_bw() + theme_boxplot

  p <- sp_ggplot_layout(
    p,
    legend.position = legend.position,
    extra_ggplot2_cmd = extra_ggplot2_cmd,
    x_label = x_label,
    y_label = y_label,
    title = title,
    additional_theme = additional_theme,
    coordinate_flip = coordinate_flip,
    ...
  )
  p

}
