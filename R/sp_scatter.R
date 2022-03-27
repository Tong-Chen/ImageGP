#' Generating scatter plot
#'
#' @param data Data file or dataframe (with header line, the first column is not the rowname, tab seperated).
#' @param xvariable The variable for x axis. NECESSARY, such `X_val` (one of column names), both text and number works.
#' @param yvariable The variable for y axis. NECESSARY, such as `Y_val` (one of column names), both text and number works.
#' @param label_variable Label points with given text. Default no-label, accept a string like
#' `Samp`` (one of column names) here to label `Samp` column text to points.
#' @param xvariable_order The order for x-axis when `xvariable`s are text. Default alphabetical order,
#' accept a string like `c('K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC')`.
#' @param yvariable_order The order for y-axis when `yvariable`s are text. Default alphabetical order,
#' accept a vector like c('K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC').
#' @param color_variable The variable for point color. Optional, such as `color` (one of column names).
#' @param color_variable_order The order for color variable. Default alphabetical order,
#' accept a vector like c('K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC').
#' @param shape_variable The variable for point shape. Optional, such as `shape` (one of column names).
#' @param shape_variable_order The order for shape variable. Default alphabetical order,
#' accept a vector like c('K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC').
#' @inheritParams sp_manual_color_ggplot2
#' @param log_variables Get log-transformed data for given variable. Default NULL, means no log10 transform.
#' Accept a vector like `c('color')` (one or several of column names) to get (-1) * log10(color).
#' @param log_transform Get log-transformed data for `log_variables`. Default `log2`, means log2 transform if
#' `log_variabels` are not null. Accept `log10`.
#' @param size_variable The variable for point size. Optional,
#' such as a number or a variable like `count` (one of column names), normally should be number column.
#' @param scale_size_min Scale size with minimum value specified
#' @param scale_size_max Scale size with maximum value specified
#' @param alpha Transparency value for points. Optional, such as a number or a variable
#' indicating one data column, normally should be number column (one of column names).
#' @param jitter Jitter points. Normally used when x and y axis variable is in text format
#' or represents group information to avoid point overlaps. Default FALSE.
#' @param jitter_text Make point labels not overlap. Default FALSE.
#' @param title Title of picture. Default empty title.
#' @param label_font_size Label font size. Default system default. Accept a number.
#' Supplying numbers less than 5 to shrink fonts.
#' @param scale_y_way The way to scale Y-axis like `scale_y_log10`, `coord_trans(y="log10")`,
#' `scale_y_continuous(trans="log2")`, `coord_trans(y="log2")`.
#' @param smooth_method The smooth method one wants to use, eg. auto, lm, glm, gam, loess, rlm.
#' For observations < 1000 default is 'loess', observations >= 1000 defaults to 'gam'.
#' Default 'no smooth' meaning show the real lines and do not smooth lines. Accept auto, lm, glm, gam, loess, rlm.
#' @param line_size line size. Default NULL. Accept a number.
#' @param facet Wrap plots by given column. This is used to put multiple plot
#' in one picture. Used when `melted` is FALSE, normally a string `set` (one of column names)
#' should be suitable for this parameter.
#' @param nrow 	The number of rows one want when `melted` is used. Default NULL.
#' @param ncol The number of columns one want when `melted` is used. Default NULL.
#' @param scales Paramter for scales for facet. Default `fixed` meaning each inner graph
#' @inheritParams sp_ggplot_layout
#' @param ... Parametes given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#' @examples
#'
#' library(plyr)
#' library(ggplot2)
#' library(grid)
#' library(ggbeeswarm)
#' library(ggrepel)
#' library(htmlwidgets)
#' library(plotly)
#' library(ImageGP)
#' scatter_test_data <- data.frame(Samp = letters[1:6], Color = sample(c("group1", "group2", "group3"),6,replace = TRUE),
#' X_val = runif(6), Y_val = runif(6), Size = sample(4:20, size = 6),
#' Shape = sample(c("cluster1","cluster2"),6,replace = TRUE))
#'
#' sp_scatterplot(data=scatter_test_data,xvariable = "X_val",yvariable = "Y_val",
#' color_variable = "Color", shape_variable= "Shape",
#' size_variable = "Size",label="Samp",Jitter = TRUE)
#'
#'
#'
#' ## Not run:
#' scatter_data = "scatter.txt"
#' sp_scatterplot(data="scatter.txt",xvariable = "X_val",yvariable = "Y_val",
#' color_variable = "Color", shape_variable= "Shape",size_variable = "Size",
#' label="Samp", xvariable_order = c(1,3,2), yvariable_order = c(2,1,3),
#' color_variable_order = c("grp2","grp1","grp3"),
#' shape_variable_order = c("cluster2","cluster1"),label_font_size=2)
#'
#' sp_scatterplot(data="scatter.txt",xvariable = "X_val",yvariable = "Y_val", color_variable = "Color", shape_variable= "Shape",
#' size_variable = "Size",label="Samp",Jitter = TRUE)
#'
#' sp_scatterplot(data="scatter.txt",xvariable = "X_val",yvariable = "Y_val", color_variable = "Color", shape_variable= "Shape",
#' size_variable = "Size",label="Samp",Jitter = TRUE,facet = "Color", scales = "free_y")
#' ## End(Not run)
#'
sp_scatterplot <- function (data,
                            xvariable = NULL,
                            yvariable = NULL,
                            label_variable = NULL,
                            xvariable_order = NULL,
                            yvariable_order = NULL,
                            color_variable_order = NULL,
                            shape_variable_order = NULL,
                            manual_color_vector = NULL,
                            log_variables = NULL,
                            log_transform = "log2",
                            size_variable = NULL,
                            geom_text_repel = TRUE,
                            shape_variable = NULL,
                            color_variable = NULL,
                            point_hjust = 0,
                            line_size = NULL,
                            smooth_method = "no smooth",
                            alpha = 1,
                            jitter = FALSE,
                            jitter_text = F,
                            scale_size_min = NULL,
                            scale_size_max = NULL,
                            coordinate_flip = FALSE,
                            legend.position = 'right',
                            xtics_angle = 0,
                            title = NULL,
                            x_label = NULL,
                            y_label = NULL,
                            extra_ggplot2_cmd = NULL,
                            label_font_size = 3,
                            facet = NULL,
                            nrow = NULL,
                            ncol = NULL,
                            scales = 'fixed',
                            scale_y_way = NULL,
                            ...) {
  if (class(data) == "character") {
    data <- sp_readTable(data)
  } else if (!"data.frame" %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  data_colnames <- colnames(data)

  if (!(xvariable %in% data_colnames &&
        yvariable %in% data_colnames)) {
    stop(paste(xvariable, 'or', yvariable, 'must be column names of data!'))
  }

  if (!is.numeric(data[[yvariable]]) &&
      !sp.is.null(yvariable_order)) {
    data = sp_set_factor_order(data, yvariable, yvariable_order)
  }

  if (!is.numeric(data[[yvariable]]) &&
      !sp.is.null(xvariable_order)) {
    data = sp_set_factor_order(data, xvariable, xvariable_order)
  }

  if (!sp.is.null(color_variable)) {
    if (!(color_variable %in% data_colnames)) {
      stop(paste(color_variable, 'must be column names of data!'))
    }

    data = sp_set_factor_order(data, color_variable, color_variable_order)

  }

  if (!sp.is.null(log_variables)) {
    for (i in log_variables) {
      y_add = sp_determine_log_add(data[[i]])
      data[[i]] <- eval(parse(text = log_transform))(data[[i]])
    }
  }

  if (!sp.is.null(shape_variable)) {
    if (!(shape_variable %in% data_colnames)) {
      stop(paste(shape_variable, 'must be column names of data!'))
    }
	  data = sp_set_factor_order(data, shape_variable, shape_variable_order)
	  data[[shape_variable]] <- as.factor(data[[shape_variable]])
	  shapes <- generate_shapes(data, shape_variable)
  }


  xvariable_en = sym(xvariable)
  yvariable_en = sym(yvariable)

  p <- ggplot(data, aes(x = !!xvariable_en, y = !!yvariable_en))

  groupOnX = NULL

  if (jitter_text) {
    geom_text = geom_text_repel
  }

  if (jitter) {
    geom_point_self <- geom_quasirandom
    groupOnX = T
    position = position_quasirandom(groupOnX = groupOnX)
  } else {
    geom_point_self <- geom_point
    position = "identity"
  }



  if (!sp.is.null(size_variable)) {
    if (numCheck(size_variable)) {
      size_variable = as.numeric(size_variable)
      p <-
        p + geom_point_self(size = size_variable,
                            alpha = alpha,
                            groupOnX = groupOnX)
    } else {
      if (!(size_variable %in% data_colnames)) {
        stop(paste(size_variable, 'must be column names of data!'))
      }
      size_variable_en = sym(size_variable)
      p <-
        p + geom_point_self(aes(size = !!size_variable_en),
                            alpha = alpha,
                            groupOnX = groupOnX)
    }
  } else {
    p <- p + geom_point_self(alpha = alpha, groupOnX = NULL)
  }

  if (!sp.is.null(color_variable)) {
    if (!(color_variable %in% data_colnames)) {
      stop(paste(color_variable, 'must be column names of data!'))
    }
    color_variable_en = sym(color_variable)
    p <- p + aes(color = !!color_variable_en)
    p <-
      sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)
  }

  if (!sp.is.null(shape_variable)) {
    shape_variable_en = sym(shape_variable)
    p <- p + aes(shape = !!shape_variable_en) +
      scale_shape_manual(values = shapes)
    }



  if (!sp.is.null(label_variable)) {
    if (!(label_variable %in% data_colnames)) {
      stop(paste(label_variable, 'must be column names of data!'))
    }
    label_variable_en = sym(label_variable)

    p <-
      p + geom_text(
        aes(label = !!label_variable_en),
        size = label_font_size,
        show.legend = F,
        position = position,
        color = "black"
      )
    if (sp.is.null(scale_size_min) && sp.is.null(scale_size_max)) {
      scale_size_min = label_font_size * 1.5
      scale_size_max = label_font_size * 4
    }
  }

  if (!sp.is.null(scale_size_min) && !sp.is.null(scale_size_max)) {
    p <- p + scale_size(name = size_variable,
                        range = range(scale_size_min, scale_size_max))
  }

  if (!sp.is.null(facet)) {
    p <- p + facet_wrap( ~  .data[[facet]],
                         nrow = nrow ,
                         ncol = ncol ,
                         scale = scales)
  }

  if (!sp.is.null(scale_y_way)) {
    p <- p + eval(parse(text = scale_y_way))
  }


  if(!sp.is.null(line_size)){
    if(numCheck(line_size)){
      line_size = as.numeric(line_size)
    } else{
      line_size_en = sym(line_size)
    }
  }

  if (smooth_method != "no smooth") {
    if(!sp.is.null(line_size) && is.numeric(line_size)){
      p <- p + stat_smooth(method = smooth_method,
                           se = FALSE,
                           size=line_size)
    }else{
      p <- p + stat_smooth(method = smooth_method,
                           se = FALSE)
    }
  }

  p <- sp_ggplot_layout(
    p,
    xtics_angle = xtics_angle,
    legend.position = legend.position,
    extra_ggplot2_cmd = extra_ggplot2_cmd,
    x_label = x_label,
    y_label = y_label,
    title = title,
    coordinate_flip = coordinate_flip,
    ...
  )
  p
}
