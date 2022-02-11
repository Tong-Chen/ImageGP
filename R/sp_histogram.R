

#' Generating histogram plot
#'
#' @param data Data frame or data file (with header line, the first
#' column will not be treated as the row-names, tab separated)
#' @param melted `TRUE` for dealing with long format matrix, the program will skip
#' melt preprocess. If input is wide format matrix, this parameter should be set to `FALSE`.
#' @param xvariable Specify the column name used as Y-axis value (one of column names,
#' should be specified when inputting long format matrix).
#' @param color_variable Specify the column name used as color variable.
#' @param color_variable_order Levels for color variable to set their appearing order.
#' Accept a vector like `c('ctcf','h3k27ac')`
#' (normally whole or a subset of unique values of one column with order specified).
#' @param group_variable Specify the column name used as group variable.
#' @param group_variable_order Levels for color variable to set their appearing order.
#' Accept a vector like `c('ctcf','h3k27ac')`
#' (normally whole or a subset of unique values of one column with order specified).
#' @param yaxis_statistics Plot with density or count or frequency in y-axis.
#' Default `frequency`, accept `density`, `count`.
#' When `plot_type` is both, `frequency`` will be given here.
#' @param plot_type Plot frequency or density or hist or both frequency and histogram.
#' Default line means frequency accept `density_line`, `hist` or `both`.
#' @param value_scale Scale value. All values will be divided by supplied one. Default 1.
#' @param color_variable_order Levels for color variable to set their appearing order.
#' Accept a vector like `c('ctcf','h3k27ac')`
#' (normally whole or a subset of unique values of one column with order specified).
#' @param facet_variable_order  Order of facets (to set the appearance order of each subplot).
#' Accept a vector like `c('ctcf','h3k27ac')`
#' (normally whole or a subset of unique values of one column with order specified).
#' @param maximum_allowed The largest value allowed. Values other than this will be set
#' as this value. Default Inf.
#' @param binwidth The value for `binwidth`` (like the width of each bin) in
#' `geom_density`, and `binwidth` for `geom_histogram`.
#' @param hist_bar_position Position parameter for hist bars. Default `identity`, accept `dodge`.
#' @param fill_area Fill the area if TRUE. Default FALSE.
#' @param line_size line size. Default 1. Accept a number.
#' @param add_mean_value_vline Add mean value as vline. Default FALSE,  accept TRUE.
#' @param xtics Show the X axis.
#' @param ytics Show the Y axis.
#' @param manual_xtics_pos Manually set the positions of xtics. Default FALSE, accept a series of
#' numbers in following format "c(1,2,3,4,5)" or other R code that can generate a vector
#' to set the position of manual xtics.
#' @param manual_xtics_value Manually set the values of xtics when `manual_xtics_position` is specified.
#' Default the content of `manual_xtics_position` when `manual_xtics_position` is specified,
#'accept a series of numbers in following format "c(1,2,3,4,5)" or other R code that can
#'generate a vector to set the value of manual xtics.

#' @inheritParams sp_ggplot_facet
#' @inheritParams sp_transfer_one_column
#' @inheritParams sp_load_font
#' @inheritParams sp_ggplot_layout
#' @inheritParams sp_load_font
#' @inheritParams sp_ggplot_add_vline_hline
#' @param ... Other parameters given to sp_ggplot_layout
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' histogram_test_data <- data.frame(Type = letters[1:2], Value = runif(80))
#' sp_histogram(data = histogram_test_data, xvariable = "Value", color_variable = "Type",
#'  plot_type = "both", yaxis_statistics = "count/sum(count)")
#'
#' ## Not run:
#' input <- "histogram.data"
#' sp_histogram(data=input, xvariable = "Value", color_variable = "Type",
#'  plot_type = "both", yaxis_statistics = "count/sum(count)")
#' ## End(Not run)

sp_histogram <- function(data ,
                         melted = FALSE,
                         xvariable = NULL,
                         color_variable = NULL,
                         color_variable_order = NULL,
                         group_variable = NULL,
                         group_variable_order = NULL,
                         value_scale = 1,
                         facet_variable = NULL,
                         facet_variable_order = NULL,
                         maximum_allowed = Inf,
                         yaxis_statistics = 'count/sum(count)',
                         plot_type = 'line',
                         binwidth = NULL,
                         hist_bar_position = "identity",
                         alpha = 0.4,
                         yaxis_scale_mode = NULL,
                         y_add = 0,
                         fill_area = FALSE,
                         line_size = 1,
                         add_mean_value_vline = FALSE,
                         manual_color_vector = "Accent",
                         facet_scales = 'fixed',
                         facet_ncol = NULL,
                         facet_nrow = NULL,
                         xtics = TRUE,
                         ytics = TRUE,
                         legend.position = "right",
                         xtics_angle = 0,
                         manual_xtics_value = NULL,
                         manual_xtics_pos = NULL,
                         custom_vline_x_position = NULL,
                         custom_vline_anno = NULL,
                         x_label = NULL,
                         y_label = NULL,
                         title = '',
                         font_path = NULL,
                         debug = F,
                         extra_ggplot2_cmd = NULL,
                         filename = NULL,
                         base_font_size = 11,
                         ...) {
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  fontname = sp_load_font(font_path = font_path)

  if (melted) {
    if (sp.is.null(group_variable) || sp.is.null(xvariable)) {
      stop("For melted matrix, <group_variable> and <xvariable> should be supplied.")
    }
  } else {
    xvariable = 'value'
    group_variable = 'variable'
  }

  data <- sp_read_in_long_wide_matrix(data, "xvariable", melted)

  #print(data)

  wide_rownames <- data$wide_rownames
  wide_colnames <- data$wide_colnames
  data <- data$data
  data_colnames <- colnames(data)


  if (sp.is.null(color_variable)) {
    color_variable = group_variable
  } else {
    data = sp_set_factor_order(data, color_variable, color_variable_order)
  }

  if (!(
    xvariable %in% data_colnames &&
    group_variable %in% data_colnames
  )) {
    stop(paste(
      group_variable,
      'or',
      xvariable,

      'must be one of column names of data!'
    ))
  }

  data = sp_set_factor_order(data, group_variable, group_variable_order)

  if (!sp.is.null(yaxis_scale_mode)) {
    data <-
      sp_transfer_one_column(
        data,
        variable = xvariable,
        yaxis_scale_mode = yaxis_scale_mode,
        y_add = y_add
      )
  }

  #stat = "density"

  #if (!is.numeric(data[[xvariable]])) {
  #  stat = "count"
  #stop("Value variable column must be numeric.")
  #}

  if (value_scale != 1) {
    data[[xvariable]] <- data[[xvariable]] / value_scale
  }


  if (!sp.is.null(facet_variable)) {
    if (!(facet_variable %in% data_colnames)) {
      stop(paste(facet_variable, 'must be one of column names of data!'))
    }
    data = sp_set_factor_order(data, facet_variable, facet_variable_order)
  }

  maximum_allowed <- as.numeric(maximum_allowed)

  if (maximum_allowed != Inf) {
    data[[xvariable]][data[[xvariable]] > maximum_allowed] <-
      maximum_allowed
  }

  xvariable_en = sym(xvariable)
  color_variable_en = sym(color_variable)
  group_variable_en = sym(group_variable)

  # print(xvariable)

  p <-
    ggplot(data, aes(
      x = !!xvariable_en,
      group = !!group_variable_en
    ))

  geom_histogram_options = list()
  geom_histogram_options$binwidth = binwidth
  geom_histogram_options <-
    geom_histogram_options[!sapply(geom_histogram_options, sp.is.null)]

  if (length(geom_histogram_options) > 1) {
    geom_histogram_options_used = T
  } else {
    geom_histogram_options_used = F
  }

  if (sp.is.null(y_label)) {
    if (yaxis_statistics == "count") {
      y_label = "Count"
    } else if (yaxis_statistics == "count/sum(count)") {
      y_label = "Frequency"
    } else if (yaxis_statistics  == "density" || plot_type == "density_line") {
      y_label = "Density"
    }
  }

  if (plot_type == "hist" ||  plot_type == "both") {
    geom_histogram_1 <- function(...) {
      geom_histogram(
        aes(
          y = after_stat(eval(parse(text = yaxis_statistics))) ,
          fill = !!color_variable_en
        ),
        alpha = alpha ,
        position = hist_bar_position ,
        bins = 30,
        stat = "bin",
        ...
      )
    }

    p <- p + do.call(geom_histogram_1, geom_histogram_options)
  }


  if (plot_type == "density_line" ||
      (plot_type  == "line" &&
       yaxis_statistics  == "density") ||
      (plot_type  == "both" &&  yaxis_statistics  == "density")) {
    geom_density_1_option = list()
    geom_density_1_option$stat = "density"
    geom_density_1_option$size = line_size
    geom_density_1_option$alpha = alpha
    if (fill_area) {
      geom_density_1 <- function(...) {
        geom_density(mapping = aes(
          fill = !!color_variable_en,
          color = !!color_variable_en
        ),
        ...)
      }
    } else {
      geom_density_1 <- function(...) {
        geom_density(mapping = aes(color = !!color_variable_en),
                     ...)
      }
    }

    p <-
      p + do.call(geom_density_1,
                  c(geom_density_1_option, geom_histogram_options))

    if (sp.is.null(y_label)) {
      y_label = "Density"
    }
  }

  if ((plot_type == "line" &&
       yaxis_statistics  != "density") ||
      (plot_type == "both" &&  yaxis_statistics  != "density")) {
    geom_density_color_fill_option = list()
    geom_density_color_fill_option$stat = "bin"
    geom_density_color_fill_option$size = line_size
    geom_density_color_fill_option$alpha = alpha

    if (fill_area) {
      geom_area_poly_1 <- function(...) {
        geom_area(mapping = aes(
          y = after_stat(eval(parse(text = yaxis_statistics))),
          fill = .data[[color_variable]],
          color = .data[[color_variable]]
        ),
        ...)
      }
    } else {
      geom_area_poly_1 <- function(...) {
        geom_area(mapping = aes(y = after_stat(eval(
          parse(text = yaxis_statistics)
        )),
        color = .data[[color_variable]]),
        ...)
      }
    }

    p <- p + do.call(geom_area_poly_1,
                     c(geom_density_color_fill_option, geom_histogram_options))
  }

  p <-
    sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)

  if (fill_area) {
    p <-
      sp_manual_fill_ggplot2(p, data, color_variable, manual_color_vector)
  }


  if (add_mean_value_vline) {
    #p <- p + ggnewscale::new_scale_color()
    group_variable_vector <- unique(c(group_variable, color_variable, facet_variable))
    group_variable_vector <- group_variable_vector[!sapply(group_variable_vector, sp.is.null)]
    data$combine__grp__for__statistis_sp <- do.call(paste0, data[group_variable_vector])

    # print(data)
    cdf <- data %>% group_by(combine__grp__for__statistis_sp) %>%
      summarise(rating.mean = mean(!!xvariable_en))
    cdf <- as.data.frame(cdf)
    rownames(cdf) <- cdf$combine__grp__for__statistis_sp

    data$rating.mean = cdf[as.character(data$combine__grp__for__statistis_sp),]$rating.mean

    #colnames(cdf) <- c(color_variable, "rating.mean")

    # print(cdf)

    p <- p + geom_vline(
      data = data,
      mapping=aes(xintercept = rating.mean,
          color = !!color_variable_en),
      linetype = "dashed",
      size = 0.5,
      show.legend = F
    ) + geom_text(
      data = data,
      mapping = aes(
        x = rating.mean * 1.01,
        y = 0,
        label = prettyNum(rating.mean, digits =
                            2),
        color = !!color_variable_en
      ),
      show.legend = F,
      hjust = 1,
      angle = 90
    )
  }


  if (!sp.is.null(facet_variable)) {
    p <-
      sp_ggplot_facet(p, facet_variable, facet_ncol, facet_nrow, facet_scales)
  }

  p <- sp_ggplot_add_vline_hline(p,
                                 custom_vline_x_position = custom_vline_x_position,
                                 custom_vline_anno = custom_vline_anno)


  additional_theme = list()
  if (!xtics) {
    additional_theme$axis.text.x = element_blank()
    additional_theme$axis.ticks.x = element_blank()
  }

  if (!ytics) {
    additional_theme$axis.text.y = element_blank()
  }

  if (!sp.is.null(manual_xtics_pos)) {
    if (sp.is.null(manual_xtics_value)) {
      manual_xtics_value <- manual_xtics_pos
    }
    p <-
      p + scale_x_continuous(breaks = manual_xtics_pos, labels = manual_xtics_value)
  }

  p <- sp_ggplot_layout(
    p,
    xtics_angle = xtics_angle,
    legend.position = legend.position,
    extra_ggplot2_cmd = extra_ggplot2_cmd,
    filename = filename,
    title = title,
    x_label = x_label,
    y_label = y_label,
    additional_theme = additional_theme,
    fontname = fontname,
    base_font_size = base_font_size,
    ...
  )
  p
}
