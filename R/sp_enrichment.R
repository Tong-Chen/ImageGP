#' Generating enrichment plot
#'
#' @param data Data file (with header line, the first column is not the rowname, tab seperated).
#' @param xvariable One of column names used as the variable for horizontal axis.
#' @param yvariable One of column names used as the variable for vertical axis.
#' @param yvariable_width Default 60. Normally yvariable would be terms. This is used to set
#' the max-allowed term length. Terms with length larger than given  value will be wrap
#' to new line.
#' @param xvariable_order The order for horizontal axis. Only woroks if horizontal axis
#' variables are strings. If horizontaol axis are numbers like GeneRatio,
#' this will be treated as shape_variable order. Default alphabetical order,
#' accept a vector like c('SampA_Up','SampB_Up').
#' @param shape_variable One of column names used as the variable for defining point shapes.
#' @param shape_variable_order The order for shape variable. Default alphabetical order,
#' accept a vector like c('SampA_Up','SampB_Up').
#' @param size_variable One of column names used as the variable for defining point sizes.
#' Normally `count`. Default the variable name given to `sqrt_transform_variable`.
#' @param scale_size_min Scale size with minimum value specified
#' @param scale_size_max Scale size with maximum value specified
#' @param title Title of picture. Default empty.
#' @param x_label X-axis title of picture. Default empty.
#' @param y_label Y-axis title of picture. Default empty.
#' @param log10_transform_variable Get log-transformed data for given variable.
#' Default `NULL`, means no log10 transform.
#' Accept a variable like `p_value` to get `(-1) * log10(p_value)`.
#' @param sqrt_transform_variable Get square root-transformed data for given variable.
#' Default `NULL`, means no sqrt transform.
#' Accept a variable like `count` to get `sqrt(count)`.
#' @param color_variable One of column names used as the variable for defining point colors. Default the
#' variable name given to `log10_transform_variable`.
#' @inheritParams sp_ggplot_layout

#' @inheritParams sp_manual_color_ggplot2
#' @param ... Parameters given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#'
#'
#' ## Not run:
#' enrichment_data <- "enrichment.data"
#'
#' library(ImageGP)
#' enrichment_data <- "enrichment.data"
#' sp_enrichment(data = enrichment_data, xvariable = "GeneRatio", yvariable = "Description",
#'               log10_transform_variable = "Qvalue", sqrt_transform_variable = "Count",
#'               shape_variable = "SampleGroup")
#'
#' enrichment_data <- "goeast.enrich.txt"
#' enrichment.data <- sp_readTable(enrichment_data)
#' head(enrichment.data)
#'
#' p <- sp_enrichment(data = enrichment.data, xvariable = "log_odds_ratio",
#'                    yvariable = "Term", color_variable = "p",
#'                    log10_transform_variable="p", size_variable = "q",
#'                    sqrt_transform_variable = "q",
#'                    shape_variable = "Ontology")
#' ## End(Not run)
#'
#'
sp_enrichment <- function(data,
                          xvariable,
                          yvariable,
                          size_variable = NULL,
                          color_variable = NULL,
                          xvariable_order = NULL,
                          yvariable_order = NULL,
                          shape_variable_order = NULL,
                          title = NULL,
                          x_label = NULL,
                          y_label = NULL,
                          scale_size_max = NULL,
                          scale_size_min = NULL,
                          yvariable_width = 60,
                          manual_color_vector = c("green", "red"),
                          log10_transform_variable = NULL,
                          sqrt_transform_variable = NULL,
                          legend.position = "right",
                          xtics_angle = 0,
                          shape_variable = NULL,
                          coordinate_flip = FALSE,
                          extra_ggplot2_cmd = NULL,
                          ...) {
  if (class(data) == "character") {
    data <- sp_readTable(data, row.names = NULL)
  } else if (class(data) != "data.frame") {
    stop("Unknown input format for `data` parameter.")
  }

  if (sp.is.null(xvariable) || sp.is.null(yvariable)) {
    stop('xvariable or yvariable must be specified!')
  }

  data_colnames <- colnames(data)

  if (!(xvariable %in% data_colnames &&
        yvariable %in% data_colnames)) {
    stop(paste(xvariable, 'or', yvariable, 'must be column names of data!'))
  }

  xval_type = "string"
  if (is.numeric(data[[xvariable]]) || numCheck(data[[xvariable]])) {
    xval_type = "numeric"
    # When meets unusual numerical type like 2/3, transfer them to numeric
    if (!is.numeric(data[[xvariable]])) {
      #print(data[[xvariable]])
      data[[xvariable]] <- mixedToFloat(data[[xvariable]])
    }
  }

  if (xval_type == "string") {
    # Do not know why add this as default, just comment out 20200705
    # if (sp.is.null(shape_variable)) {
    #   shape_variable = xvariable
    # }
    data = sp_set_factor_order(data, xvariable, xvariable_order)
  }



  if (!sp.is.null(shape_variable)) {
    if (shape_variable != xvariable) {
      data = sp_set_factor_order(data, shape_variable, shape_variable_order)
    } else if (sp.is.null(xvariable_order)) {
      if (!sp.is.null(shape_variable_order)) {
        data = sp_set_factor_order(data, shape_variable, shape_variable_order)
      }
    }

    shape_level <- length(unique(data[[shape_variable]]))
    shapes = (1:shape_level) %% 30

    shape_variable_en = sym(shape_variable)
  }

  # First order by Term, then order by shape_variable
  if (!sp.is.null(shape_variable) & xval_type != "numeric") {
    data <- data[order(data[[yvariable]], data[[shape_variable]]),]
  }

  if (!sp.is.null(color_variable) &&
	  !is.numeric(data[[color_variable]]) &&
      numCheck(data[[color_variable]])) {
    data[[color_variable]] = mixedToFloat(data[[color_variable]])
  }

  if (!sp.is.null(log10_transform_variable) &&
      (sp.is.null(color_variable) ||
       color_variable == log10_transform_variable)) {
    color_variable = paste0("negLog10_", log10_transform_variable)
  }

  if (!sp.is.null(log10_transform_variable)) {
    log_name = paste0("negLog10_", log10_transform_variable)
    col_name_data <- colnames(data)
    col_name_data <- c(col_name_data, log_name)
    if (!numCheck(data[[log10_transform_variable]])) {
      stop(
        paste(
          log10_transform_variable,
          "column is not numerical column. Please do not set log10 transform on this column."
        )
      )
    } else {
      if (!is.numeric(data[[log10_transform_variable]])) {
        data[[log10_transform_variable]] = mixedToFloat(data[[log10_transform_variable]])
      }
    }
    data$log_name <- log10(data[[log10_transform_variable]]) * (-1)
    data$log_name[data$log_name == Inf] = max(data$log_name[data$log_name !=
                                                              Inf]) + 2
    colnames(data) <- col_name_data
  }

  if (!sp.is.null(size_variable) &&
      numCheck(data[[size_variable]]) &&
      !is.numeric(data[[size_variable]])) {
    data[[size_variable]] = mixedToFloat(data[[size_variable]])
  }

  if (!sp.is.null(sqrt_transform_variable) &&
      (sp.is.null(size_variable) ||
       size_variable == sqrt_transform_variable)) {
    size_variable = paste0("sqrt_", sqrt_transform_variable)
  }

  if (!sp.is.null(sqrt_transform_variable)) {
    sqrt_name = paste0("sqrt_", sqrt_transform_variable)
    col_name_data <- colnames(data)
    col_name_data <- c(col_name_data, sqrt_name)
    if (!numCheck(data[[sqrt_transform_variable]])) {
      stop(
        paste(
          sqrt_transform_variable,
          "column is not numerical column. Plase do not set sqrt transform on this column."
        )
      )
    } else {
      if (!is.numeric(data[[sqrt_transform_variable]])) {
        data[[sqrt_transform_variable]] = mixedToFloat(data[[sqrt_transform_variable]])
      }
    }
    data$sqrt_name <- sqrt(data[[sqrt_transform_variable]])

    colnames(data) <- col_name_data
  }

  if (sp.is.null(yvariable_order)) {
    # Get the count of each unique Term
    data_freq <- as.data.frame(table(data[[yvariable]]))

    colnames(data_freq) <- c(yvariable, "IDctct")

    data2 <- merge(data, data_freq, by = yvariable)

    if (!sp.is.null(shape_variable)) {
      # Collapse shape_variable for each Term

      data_samp <-
        aggregate(as.formula(paste(shape_variable, "~", yvariable)),
                  data2, paste, collapse = "_")
      colnames(data_samp) <- c(yvariable, "sam_ct_ct_ct")
      data2 <- merge(data2, data_samp, by = yvariable)

      #print(data2)

      if (xval_type != "string") {
        data3 <-
          data2[order(data2$IDctct,
                      data2$sam_ct_ct_ct,
                      data2[[shape_variable]],
                      data2[[xvariable]],
                      data2[[color_variable]]),]
      } else {
        data3 <-
          data2[order(data2$IDctct, data2$sam_ct_ct_ct, data2[[shape_variable]],
                      data2[[color_variable]]),]
      }
    } else{
      if (xval_type != "string") {
        data3 <-
          data2[order(data2$IDctct, data2[[xvariable]], data2[[color_variable]]),]
      } else {
        data3 <- data2[order(data2$IDctct, data2[[color_variable]]),]
      }
    }
    #print(data3)

    yvariable_order <- unique(data3[[yvariable]])
    rm(data_freq, data2, data3)
  }

  data = sp_set_factor_order(data, yvariable, yvariable_order)

  xvariable_en = sym(xvariable)
  yvariable_en = sym(yvariable)
  size_variable_en = sym(size_variable)
  color_variable_en = sym(color_variable)

  p <-
    ggplot(data, aes(x = !!xvariable_en, y = !!yvariable_en)) + geom_point()

  if (!sp.is.null(shape_variable)) {
    p <- p + aes(shape = !!shape_variable_en)
    if (shape_level > 6) {
      p <- p + scale_shape_manual(values = shapes)
    }
  }

  if (!sp.is.null(size_variable)) {
    p <- p + aes(size = !!size_variable_en)
    if(is.numeric(data[[size_variable]])){
      if(all(data[[size_variable]] == as.integer(data[[size_variable]]))){
        min = min(data[[size_variable]])
        max = max(data[[size_variable]])

		# 4 is length
		step = ceiling((max-min)/4)

        p <- p + scale_size_continuous(breaks=seq(min, max, by=step))
      }
    }
    if (!sp.is.null(scale_size_min) && !sp.is.null(scale_size_max)) {
      p <- p + scale_size(name = size_variable,
                          range = range(scale_size_min, scale_size_max))
    }
  }

  if (!sp.is.null(color_variable)) {
    p <- p + aes(color = !!color_variable_en)
    p <-
      sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)
  }


  p <- p + scale_y_discrete(
    labels = function(x)
      str_wrap(x, width = yvariable_width)
  )



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
