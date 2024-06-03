#' Generating lines plot
#'
#' @param data Data frame or data file (with header line, the first column will not be treated as
#' the rowname, tab seperated)
#' @param melted When true, it will skip melt preprocesses and each column wolud be treated as separate attributes.
#' Default FALSE, accept TRUE.
#' @param xvariable Name for x-axis variable. For normal matrix, default the first column will be used,
#' The program will assign an value 'xvariable' to represent it.
#' @param xvariable_order Set orders of x-variable. For normal matrix, `xvariable_order` is the row order of data,
#' For melted matrix, accept a vector like c('CTCF','H3K27ac','Enhancer') to set your own order.
#' @param yvariable Name for y-axis variable. For normal matrix, all numbers except the first row and first column
#' would be used. The program will assign an value 'value' to represent it. Default value.
#' @param y_add A number to add if scale is used. Default 0 meaning the minimum non-zero value wolud be used.
#' @param yaxis_scale_mode Give the following `scale_y_log10()`,
#' `coord_trans(y="log10")`, or other legal command for ggplot2 or simply `log2` to set the scale way.
#' @param legend_variable Name for legend variable. Default variable for normal matrix (un-melted),
#' When `melted` is TRUE, this parameter is required unless all points belong to one group.
#' @param legend_variable_order Set orders of legend variable. Default column order for normal matrix,
#' accept a vector like c('CTCF','H3K27ac','Enhancer') to set your own order.
#' @param color_variable Name for color variable. Normally used when wanting to color lines by groups.
#' Default same as `variable`, this should only be set when `melted` is TRUE.
#' @param color_variable_order Set orders of color variable.
#' Only used when `color_variable` is not the same as `legend_variable`.
#' @param y_start_from_zero Set Y-axis ranges from zero. Default TRUE.
#' @param x_level Levels for x-axis variable, suitable when x-axis is not treated as numerical.
#' Default the order of first column for normal matrix.
#' @param x_label Xlab label. Default empty.
#' @param y_label Ylab label. Default empty.
#' @param smooth_method The smooth method one wants to use, eg. auto, lm, glm, gam, loess, rlm.
#' For observations < 1000 default is 'loess', observations >= 1000 defaults to 'gam'.
#' Default 'no smooth' meaning show the real lines and do not smooth lines. Accept auto, lm, glm, gam, loess, rlm.
#' @param line_size line size. Default NULL. Accept a number.
#' @param alpha Color transparency value. Default 0.6. Accept a number between 0-1.
#' @inheritParams sp_manual_color_ggplot2
#' @param xtics Default TRUE to display xtics. Accept FALSE to turn off xtics.
#' @param xtics_angle Rotation angle for x-tics (anti clockwise), Default 0.
#' @param ytics Default TRUE to display ytics. Accept FALSE to turn off ytics.
#' @param manual_xtics_pos Manually set the positions of xtics. Default FALSE,
#' accept a vector of numbers in following format "c(1,2,3,4,5)" or other R code
#' that can generate a vector to set the position of xtics.
#' @param manual_xtics_value Manually set the values of xtics when `xtics_pos` is specified.
#' Default the content of `xtics_pos` when `xtics_pos` is specified,
#' accept a vector of numbers in following format "c(1,2,3,4,5)" or other R code
#' that can generate a vector to set the values of xtics.
#' @param filename Output name of pictures.
#' @param facet_wrap_formula The formula for `facet_wrap`. Default NULL meaning no facets. Currently unused.
#' @param facet_wrap_ncol Number of facet columns. Default NULL meaning no facets. Currently unused.
#' @param facet_wrap_nrow Number of facet rows. Default NULL meaning no facets. Currently unused.
#' @inheritParams sp_ggplot_add_vline_hline
#' @inheritParams sp_ggplot_layout
#' @inheritParams sp_barplot
#' @param ... Parametres given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' res_output <- data.frame(Pos=1:10,value =runif(20))
#' value=0.5
#' res_output$Variable <- ifelse(res_output$value<=value,"groupA", "groupB")
#' sp_lines(data=res_output, color_variable = "Variable", variable="Variable")
#'
#'
#' ## Not run:
#' lines_data = "lines.data"
#'
#' sp_lines(data=lines_data, xvariable = "Pos", yvariable = "value", color_variable = "Variable", variable="Variable")
#' ## End(Not run)
#'
sp_lines <- function(data,
                     melted = FALSE,
                     xvariable = NULL,
                     xvariable_order = NULL,
                     y_add = 0,
                     yvariable = NULL,
                     yaxis_scale_mode = "",
                     legend_variable = NULL,
                     legend_variable_order = NULL,
                     color_variable = NULL,
                     color_variable_order = NULL,
                     error_bar_variable = NULL,
                     x_label = NULL,
                     y_label = NULL,
                     y_start_from_zero = T,
                     title = "",
                     smooth_method = "no smooth",
                     line_size = NULL,
                     alpha = 0.6,
                     manual_color_vector = NULL,
                     xtics = TRUE,
                     xtics_angle = 0,
                     ytics = TRUE,
                     legend.position = "right",
                     manual_xtics_pos = NULL,
                     manual_xtics_value = NULL,
                     width = 10,
                     height = 10,
                     custom_vline_x_position = NULL,
                     custom_vline_anno = NULL,
                     custom_hline_y_position = NULL,
                     custom_hline_anno = NULL,
                     facet_wrap_formula = NULL,
                     facet_wrap_nrow = NULL,
                     facet_wrap_ncol = NULL,
                     filename = NULL,
                     coordinate_flip = FALSE,
                     debug = FALSE,
                     ...) {

  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if (melted) {
    if (sp.is.null(xvariable) || sp.is.null(yvariable)) {
      stop("For melted matrix, <xvariable> and <yvariable> should be supplied.")
    }
  } else {
    xvariable = 'xvariable'
    yvariable = 'value'
    legend_variable = 'variable'
  }

  if(sp.is.null(xvariable) || sp.is.null(yvariable)){
    stop('xvariable or yvariable must be specified!')
  }

  data <- sp_read_in_long_wide_matrix(data, xvariable, melted)

  #print(data)

  wide_rownames <- data$wide_rownames
  wide_colnames <- data$wide_colnames
  data <- data$data
  data_colnames <- colnames(data)

  if (!melted){
    xvariable_order = wide_rownames
    color_variable_order = wide_colnames
  }

    if (!(xvariable %in% data_colnames &&
          yvariable %in% data_colnames)) {
      stop(paste(xvariable, 'or', yvariable, 'must be column names of data!'))
    }


  if (sp.is.null(legend_variable)) {
    cat("All points would be treated as in one group.\n")
    legend_variable = '_all_in_one_grp_sp_'
    data[[legend_variable]] = 'Group'
  }

  if (yaxis_scale_mode != "") {
    # Give the minimum non-zero value to add to avoid log2(0)
    if (y_add == 0) {
      y_add = sp_determine_log_add(data[[yvariable]])
    }
    data[[yvariable]] <- data[[yvariable]] + y_add
    if (yaxis_scale_mode == "log2") {
      data[[yvariable]] <- log2(data[[yvariable]])
    }
    if (yaxis_scale_mode == "log10") {
      data[[yvariable]] <- log10(data[[yvariable]])
    }
  }


  xval_type = "string"
  if (numCheck(data[[xvariable]])) {
    xval_type = "numeric"
  }

  data = sp_set_factor_order(data, xvariable, xvariable_order)


  if (sp.is.null(color_variable)) {
    color_variable = legend_variable
  }

  if (color_variable != "variable" &&
      color_variable != legend_variable) {
    data = sp_set_factor_order(data, color_variable, color_variable_order)
  }


  if (!sp.is.null(legend_variable_order)) {
    data = sp_set_factor_order(data, legend_variable, legend_variable_order)
  } else if (!melted) {
    data[[legend_variable]] <-
      factor(data[[legend_variable]], levels = wide_colnames, ordered = TRUE)
  }

    xvariable_en = sym(xvariable)
  yvariable_en = sym(yvariable)

  color_variable_en = sym(color_variable)
  legend_variable_en = sym(legend_variable)

  p <-
    ggplot(
      data,
      aes(
        x = !!xvariable_en,
        y = !!yvariable_en,
        color = !!color_variable_en,
        group = !!legend_variable_en
      )
    )
  # +geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1, size=0.5)

  if (!sp.is.null(error_bar_variable)) {
    if (!(error_bar_variable %in% c(data_colnames, "Standard_deviation"))) {
      stop(paste(error_bar_variable,'must be column names of data!'))
    }


  # legend_variable_en = sym(legend_variable)
  legend_variable_vector <- unique(c(xvariable, legend_variable))
  legend_variable_vector <- legend_variable_vector[!sapply(legend_variable_vector, sp.is.null)]
  if (length(legend_variable_vector) == 1 ){
    xvariable = legend_variable_vector
    color_variable = legend_variable_vector
  } else {
    color_variable = legend_variable
  }
  # library(Rmisc)
  # summarySE(data, measurevar="Binding_strength", groupvars=c("Epi_type","Pos"))
  # data$Pos <- as.character(data$Pos)
  data_sd_mean <- data %>% group_by(across(legend_variable_vector)) %>%
    summarise(Standard_deviation=sd(!!yvariable_en), Mean_value=mean(!!yvariable_en)) %>%
    ungroup() %>%
    group_by(!!xvariable_en) %>%
    mutate(Mean_value_cumsum_s_p=rev(cumsum(rev(Mean_value))))

  data <- as.data.frame(data_sd_mean)
  print(data_sd_mean)

  # data_sd_mean = sp_set_factor_order(data_sd_mean, legend_variable, legend_variable_order)


  # data <- merge(data, data_sd_mean, by=legend_variable, all=F)

  yvariable = "Mean_value"
  yvariable_en = sym(yvariable)

  error_bar_variable = "Standard_deviation"
  error_bar_variable_en = sym(error_bar_variable)

  if(sp.is.null(color_variable)){
    color_variable <- legend_variable[legend_variable!=xvariable][1]
  }
    errorbar_base_variable = yvariable
    errorbar_base_variable_en = sym(errorbar_base_variable)
    if(!sp.is.null(group_variable)){
      p <-
        p + geom_errorbar(
          mapping = aes(
            ymin = !!errorbar_base_variable_en - !!error_bar_variable_en,
            ymax = !!errorbar_base_variable_en + !!error_bar_variable_en,
            group=!!color_variable_en
          ),
          data = data_sd_mean,
          colour = "black",
          width = 0.2,
          position = "identity"
          #position = position
        )
    } else {
      p <-
        p + geom_errorbar(
          aes(
            ymin = !!errorbar_base_variable_en - !!error_bar_variable_en,
            ymax = !!errorbar_base_variable_en + !!error_bar_variable_en,
            group=!!color_variable_en
          ),
          colour = "black",
          width = 0.2,
          position = "identity"
          #position = position
        )
    }
  }


  if (y_start_from_zero) {
    p <- p + expand_limits(y = 0)
  }

  #p <- p + scale_y_continuous(expand=c(0, 0))

  p <- p + theme(legend.key = element_blank())

  # auto compute width and height/

  if (width == 0 || height == 0) {
    if (xval_type == "string") {
      x_len = length(unique(data[[xvariable]]))
      if (is.na(x_len)) {
        x_len = 5
      }
      total_len = sum(nchar(unique(as.character(data[[xvariable]]))))
      if (is.na(total_len)) {
        total_len = 15
      }
      average_len = total_len / x_len
      if (is.na(average_len)) {
        average_len = 9
      }
      if (average_len < 10) {
        width = total_len / 3
        if (width < x_len) {
          width = 10
        }
      } else {
        if (x_len < 10) {
          width = 10
        } else if (x_len < 20) {
          width = 12 + (x_len - 8) / 5
        } else if (x_len < 100) {
          width = 14 + (x_len - 20) / 5
          if (width > 30) {
            width = 30
          }
        } else {
          width = 30
        }
      }
    } else {
      # span = max(data[[xvariable]]) - min(data[[xvariable]])
      # width = span / 1000
      width = 10
      #if (width < 10) {
      #  width = 10
      #}
    }

    if (width < 6) {
      width = 6
    }
    #height = 0.75 * width
    if (legend.position %in% c("left", "right")) {
      legend_len = max(sapply(as.vector(unique(data[[color_variable]])), nchar))
      #print(legend_len)
      if (is.na(legend_len)) {
        legend_len = 6
      }
      #if (average_len > 1) {
      height = legend_len / 8 * width + average_len / 5
      #print(height)
      if (height < 0.75 * width) {
        height = 0.75 * width
      }
      width = width + legend_len / 2
      #} else {
      #	height = legend_len / 8 * width
      #}
    } else if (legend.position %in% c("top", "bottom")) {
      total_len = sum(nchar(unique(as.character(data[[color_variable]]))))
      if (total_len / width > 4) {
        width = width * total_len / width / 4
      }
      height = width
    }
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
                           alpha = alpha, linewidth=line_size)
    }else{
      p <- p + stat_smooth(method = smooth_method,
                           se = FALSE,
                           alpha = alpha)
    }

  } else{
    if(!sp.is.null(line_size) && is.numeric(line_size)){
      p <- p + geom_line(alpha = alpha, linewidth=line_size)
    }
    else {
      p <- p + geom_line(alpha = alpha)
    }
  }


  if (!sp.is.null(line_size) && !is.numeric(line_size)) {
      # cat("Generate line size column.\n")
      if (!line_size %in% colnames(data)) {
        stop("Unexisted columns specifed for line_size parameters")
      }
      p <-
        p + aes(linewidth = !!line_size_en)

  }

  if (yaxis_scale_mode != "" && (yaxis_scale_mode  != "log2") && (yaxis_scale_mode  != "log10")) {
    # Transfer string to R code
    p <- p + eval(parse(text=yaxis_scale_mode))
  }

  p <- sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)


  # p <- p + facet

  p <- sp_ggplot_add_vline_hline(
    p,
    custom_vline_x_position = custom_vline_x_position,
    custom_hline_y_position = custom_hline_y_position,
    custom_vline_anno = custom_vline_anno,
    custom_hline_anno = custom_hline_anno
  )


  if (!xtics) {
    p <-
      p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  if (!ytics) {
    p <-
      p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  if (!sp.is.null(manual_xtics_pos)) {
    if (sp.is.null(manual_xtics_value)) {
      manual_xtics_value <- manual_xtics_pos
    }
    p <-
      p + scale_x_continuous(breaks = manual_xtics_pos, labels = manual_xtics_value)
  }

  p <- sp_ggplot_layout(p,
                        filename = filename,
                        xtics_angle = xtics_angle,
                        legend.position = legend.position,
                        x_label = x_label,
                        y_label = y_label,
                        title = title,
                        coordinate_flip = coordinate_flip,
                        width = width,
                        height = height,
                        ...)
  p
}
