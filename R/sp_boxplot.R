#' Generating box plot
#'
#' `metadata`
#'
#' @param data Data file (with header line, the first row is the colname,
#' tab separated. Multiple formats are allowed and described above)
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
#' @param legend_variable The column represents the legend information.
#' Default `xvariable` if not specified.
#' @param group_variable_for_line Specify the group of points to line together (one column name).
#' @param group_variable_order_for_line Levels for group variable for lines.
#' @param legend_variable_order Levels for legend variable.
#' Default data order, accept a vector like c('TP16','TP22','TP23') for `legend_variable` column.
#' @param legend_variable_cut Self-define intervals for legend variable when
#' values in `legend_variable` column is continuous numbers.
#' @param xvariable_order xvariable_order Levels for x-axis variable. Default data order, accept input like c('g','a','j','x','s','c','o','u') for Set column.
#' @param xvariable_cut sxvariable_cut Self-define intervals for x-axis variable.
#' @inheritParams sp_transfer_one_column
#' @param notch Using notch (sand clock shape) or not. Default FALSE.
#' @param outlier Exclude outliers. Exclude outliers or not, default `FALSE`` means keeping outliers.
#' @param out_scale The scales for one want to set to exclude outliers.
#' Default 1.05. No recommend to change unless you know what you are doing.
#' @param violin Do violin plot plus inner boxplot.
#' @param violin_nb Do violin plot without inner boxplot.
#' @param scale_violin The value given to scale for violin plot.
#' if "area", all violins have the same area (before trimming the tails).
#' If "count", areas are scaled proportionally to the number of observations.
#' If "width", all violins have the same maximum width. 'equal' is also accepted.
#' @param ID_var Other columns one want to treat as ID variable columns
#' except the one given to `xvariable`.
#' @param jitter Do jitter plot instead of boxplot.
#' @param jitter_bp Do jitter plot overlay with violin plot or boxplot or both.
#' @param dotplot Do dotplot plot instead of boxplot.
#' @param dotplot_bp Do dotplot plot overlay with violin plot or boxplot or both.
#' @param coordinate_flip Rotate the plot from vertical to horizontal.
#' Usefull for plots with many values or very long labels at X-axis
#' @param facet_variable_order The levels of wrapping to set the order of each group.
#' @param facet_singlecell_style Use specified style for showing single cell gene expression profile. Default FALSE.
#' @inheritParams sp_ggplot_facet
#' @inheritParams sp_ggplot_layout
#' @inheritParams sp_manual_color_ggplot2
#' @param ... Parametes given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' box_test_data <- data.frame(ID = letters[1:4],
#' Gene = letters[c(8,8,9,9,10,10,11,11)], Expr = runif(16))
#' sp_boxplot(data = box_test_data, xvariable = "ID", value = "Expr", variable ="Gene")
#'
#' ## Not run:
#' box_data = "box.data"
#'
#' sp_boxplot(data = box_data, xvariable = "Gene", value = "Expr", variable="Group")
#' ## End(Not run)
#'
sp_boxplot <- function(data,
                       melted = FALSE ,
                       xvariable = NULL,
                       yvariable = NULL,
                       legend_variable = NULL,
                       statistics = FALSE,
                       xtics_angle = 0,
                       legend_variable_order = NULL,
                       legend_variable_cut = NULL,
                       xvariable_order = NULL,
                       xvariable_cut = NULL,
                       group_variable_for_line = NULL,
                       group_variable_order_for_line = NULL,
                       y_add = 0,
                       yaxis_scale_mode = NULL,
                       notch = FALSE,
                       par = NULL,
                       outlier = FALSE,
                       out_scale = 1.05,
                       legend.position = 'right',
                       manual_color_vector = NULL,
                       violin = FALSE,
                       violin_nb = FALSE,
                       scale_violin = 'width',
                       ID_var = c(),
                       jitter = FALSE ,
                       jitter_bp =  FALSE ,
                       dotplot = FALSE,
                       dotplot_bp = FALSE,
                       colormodel = 'srgb',
                       coordinate_flip = FALSE,
                       facet_variable = NULL,
                       facet_variable_order = NULL,
                       x_label = NULL,
                       y_label = NULL,
                       title = NULL,
                       facet_nrow = NULL,
                       facet_ncol = NULL,
                       facet_singlecell_style = F,
                       facet_scales = 'fixed',
                       metadata = NULL,
                       debug = F,
                       filename = NULL,
                       extra_ggplot2_cmd = NULL,
                       ...) {
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if (!melted) {
    if (sp.is.null(yvariable)) {
      yvariable = "value"
    }
    if (sp.is.null(legend_variable)) {
      legend_variable = "variable"
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
  } else if ("data.frame" %in% class(data)) {
    if (!melted) {
      if (sp.is.null(xvariable)) {
        xvariable = colnames(data)[1]
      }
    }
  } else if (!"data.frame" %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  if (sp.is.null(xvariable) || sp.is.null(yvariable)) {
    stop('xvariable or yvariable must be specified!')
  }

  if (sp.is.null(legend_variable)) {
    legend_variable = xvariable
  }

  # print(data)
  if (!sp.is.null(metadata)) {
    if (class(metadata) == "character") {
      metadata <- sp_readTable(metadata, row.names = NULL)
    } else if (!"data.frame" %in% class(data)) {
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

    data[[matched_column[2]]] = data[[matched_column[1]]]
  }

  data_colnames <- colnames(data)

  if (!(xvariable %in% data_colnames &&
        yvariable %in% data_colnames)) {
    stop(paste(xvariable, 'or', yvariable, 'must be column names of data!'))
  }



  if (!sp.is.null(yaxis_scale_mode)) {
    data <-
      sp_transfer_one_column(
        data,
        variable = yvariable,
        yaxis_scale_mode = yaxis_scale_mode,
        y_add = y_add
      )
  }

  if (!sp.is.null(legend_variable_cut)) {
    data[[legend_variable]] <-
      cut(data[[legend_variable]], legend_variable_cut)
  } else if (!sp.is.null(legend_variable_order)) {
    data = sp_set_factor_order(data, legend_variable, legend_variable_order)
  }

  if (!sp.is.null(group_variable_for_line) &&
      !sp.is.null(group_variable_order_for_line)) {
    data = sp_set_factor_order(data,
                               group_variable_for_line,
                               group_variable_order_for_line)
  }
  # May be not needed. Comment out first.
  # data[[legend_variable]] <- as.factor(data[[legend_variable]])

  # A little more complex logic
  # xvariable != legend_variable
  # xvariable == legend_variable and legend not cut
  if (!sp.is.null(xvariable_cut) &&
      ((xvariable != legend_variable) ||
       sp.is.null(legend_variable_cut))) {
    data[[xvariable]] <- cut(data[[xvariable]], xvariable_cut)
  } else if (!sp.is.null(xvariable_order)) {
    data = sp_set_factor_order(data, xvariable, xvariable_order)
  }

  data[[xvariable]] <- as.factor(data[[xvariable]])

  if (!sp.is.null(facet_variable)) {
    data = sp_set_factor_order(data, facet_variable, facet_variable_order)
  }

  xvariable_en = sym(xvariable)
  yvariable_en = sym(yvariable)
  legend_variable_en = sym(legend_variable)

  p <- ggplot(data, aes(!!xvariable_en, !!yvariable_en))

  outlier.colour = 'red'



  if (!sp.is.null(group_variable_for_line) || jitter || jitter_bp) {
    outlier.colour = 'NA'
  }

  if (violin) {
    p <- p + geom_violin(
      aes(fill = !!legend_variable_en),
      stat = "ydensity",
      trim = TRUE,
      scale = scale_violin,
      position = position_dodge(width = 0.9)
    ) +
      geom_boxplot(
        aes(fill = !!legend_variable_en),
        alpha = .25,
        width = 0.15,
        position = position_dodge(width = .9),
        outlier.colour = 'NA'
      ) +
      stat_summary(
        aes(group = !!legend_variable_en),
        fun = mean,
        geom = "point",
        fill = "black",
        shape = 19,
        size = 1,
        position = position_dodge(width = .9)
      )
  } else if (violin_nb) {
    p <- p + geom_violin(
      aes(fill = !!legend_variable_en),
      stat = "ydensity",
      #position = "dodge",
      trim = TRUE,
      scale = scale_violin,
      position = position_dodge(width = .9)
    )
  } else if (jitter) {
    p <-
      p + geom_quasirandom(
        aes(
          colour = !!legend_variable_en,
          group = !!legend_variable_en
        ),
        groupOnX = NULL,
        position = position_dodge(width = .9)
      )
    p <-
      p + stat_summary(
        fun = mean,
        geom = "text",
        label = "+",
        size = 10,
        color = "black",
        position = position_dodge(width = .9)
      )
  } else if (dotplot) {
    p <- p + geom_dotplot(
      binaxis = 'y',
      aes(
        group = !!legend_variable_en,
        fill = !!legend_variable_en
      ),
      position = position_dodge(width = .9),
      stackdir = 'center',
      stackratio = 1.5,
      binwidth = .1,
      binpositions = "all",
      dotsize = 0.6,
      alpha = .75,
      na.rm = TRUE
    )
  } else {
    if (notch) {
      if (outlier) {
        p <- p + geom_boxplot(
          aes(fill = !!legend_variable_en),
          notch = TRUE,
          notchwidth = 0.3,
          outlier.colour = 'NA'
        )
      } else{
        p <- p + geom_boxplot(
          aes(fill = !!legend_variable_en),
          notch = TRUE,
          outlier.colour = outlier.colour,
          notchwidth = 0.3,
          position = position_dodge(width = .9)
        )
      }
    } else {
      if (outlier) {
        p <- p + geom_boxplot(
          aes(fill = !!legend_variable_en),
          outlier.colour = 'NA',
          position = position_dodge(width = .9)
        )
      } else{
        p <- p + geom_boxplot(
          aes(fill = !!legend_variable_en),
          outlier.colour = outlier.colour,
          position = position_dodge(width = 0.9)
        )
      }
    }
  }

  if (jitter_bp) {
    p <-
      p + geom_quasirandom(
        # p + geom_point(
        aes(group = !!legend_variable_en),
        varwidth = T,
        groupOnX = TRUE,
        dodge.width = 0.9,
        position = position_dodge(width = .9)
      )
  }
  if (dotplot_bp) {
    p <- p + geom_dotplot(
      binaxis = 'y',
      aes(group = !!legend_variable_en),
      position = position_dodge(width = .9),
      stackdir = 'center',
      stackratio = 1.5,
      binwidth = .1,
      binpositions = "all",
      dotsize = 0.6,
      alpha = .75,
      fill = "lightseagreen",
      colour = "lightseagreen",
      na.rm = TRUE
    )
  }

  if (!sp.is.null(group_variable_for_line)) {
    if (!(group_variable_for_line %in% data_colnames)) {
      stop(paste(group_variable_for_line, 'must be column names of data!'))
    }
    group_variable_for_line_en = sym(group_variable_for_line)

    if (!(jitter || jitter_bp)) {
      p <-
        p + geom_quasirandom(
          aes(group = !!legend_variable_en),
          varwidth = T,
          groupOnX = TRUE,
          dodge.width = 0.9,
          position = position_dodge(width = .9)
        )
    }

    p <- p + geom_line(
      aes(
        group = !!group_variable_for_line_en,
        color = !!group_variable_for_line_en
      ),
      position = position_quasirandom()
    )
  }


  if (!sp.is.null(yaxis_scale_mode) &&
      (yaxis_scale_mode  != "log2") &&
      (yaxis_scale_mode  != "log10")) {
    p <- p +  eval(parse(text = yaxis_scale_mode))
    # Do not know why add this
    # p <-
    #   p + stat_summary(
    #     fun.y = "mean",
    #     geom = "text",
    #     label = "----",
    #     size = 10,
    #     color = "black"
    #   )
  }

  p <-
    sp_manual_fill_ggplot2(p, data, legend_variable, manual_color_vector)

  if (statistics) {
    if (!sp.is.null(group_variable_for_line)) {
      p <- p + ggnewscale::new_scale_color()
    }
    # 代码修改自 amplicon包 microbiota/amplicon
    # https://github.com/microbiota/amplicon/blob/master/R/alpha_boxplot.R

    # data$combine_xvariable <- paste0(data[[variable]],"_",data[[xvariable]])
    # model = aov(data[[value]] ~ data$combine_xvariable, data = data)
    # if (length(unique(data[[xvariable]])) == 2) {
    #   library(agricolae)
    #   out = LSD.test(model, "data[[xvariable]]", p.adj = "none")
    #   stat = out$groups
    #   data$stat = stat[as.character(data$combine_xvariable), ]$groups
    # } else{
    #   Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
    #   Tukey_HSD_table = as.data.frame(Tukey_HSD$`data$combine_xvariable`)
    #   Tukey.levels = Tukey_HSD$`data$combine_xvariable`[, 4]
    #   Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
    #   Tukey.labels$group = rownames(Tukey.labels)
    #   Tukey.labels = Tukey.labels[order(Tukey.labels$group),]
    #   data$stat = Tukey.labels[as.character(data$combine_xvariable), ]$Letters
    # }
    #
    # max=max(data[,c(value)])
    # min=min(data[,value])
    # x = data[,c('combine_xvariable',value)]
    # y = x %>% group_by(combine_xvariable) %>% summarise_(Max=paste('max(',value,')',sep=""))
    # y=as.data.frame(y)
    # colnames(y) <- c("group","Max")
    # rownames(y)=y$group
    # data$y=y[as.character(data$combine_xvariable),]$Max + (max-min)*0.05
    #
    # x1=factor(data$combine_xvariable,levels=c(unique(data$combine_xvariable)))
    # p <- ggplot(data, aes(x1,!!value_en))
    #
    # p <- p + geom_text(data=data, aes(x=combine_xvariable, y=y,
    #                              color=!!xvariable_en, label=stat)) +
    #   theme(axis.text.x=element_text(angle =90,vjust=0.3))

    # if (xvariable != legend_variable) {
    #   data$combine__grp__for__statistis_sp <-
    #     paste(data[[xvariable]], data[[legend_variable]], sep = "___")
    # } else {
    #   data$combine__grp__for__statistis_sp <- data[[xvariable]]
    # }
    # if (!sp.is.null(facet_variable)) {
    #   data$combine__grp__for__statistis_sp <-
    #     paste(data$combine__grp__for__statistis_sp, data[[facet_variable]], sep = "___")
    # }

    group_variable_vector <- unique(c(xvariable, legend_variable, facet_variable))
    group_variable_vector <- group_variable_vector[!sapply(group_variable_vector, sp.is.null)]
    #data2 <- data[,group_variable_vector]
    data$combine__grp__for__statistis_sp <- do.call(paste0, data[group_variable_vector])

    formula = as.formula(paste(yvariable, "~", "combine__grp__for__statistis_sp"))
    # model = aov(data[[yvariable]] ~ data[[xvariable]], data = data)
    # print(formula)
    model = aov(formula, data = data)
    # print(model)
	Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
	# return(Tukey_HSD)
	Tukey_HSD_table = as.data.frame(Tukey_HSD$combine__grp__for__statistis_sp)
	suppressWarnings(sp_writeTable(Tukey_HSD_table, file = paste0(filename, ".significance.txt")))

    if (length(unique(data$combine__grp__for__statistis_sp)) == 2) {
       library(agricolae)
       out = LSD.test(model, "combine__grp__for__statistis_sp", p.adj = "none")
       # print(out)
       LSD.test_table = as.data.frame(out$statistics)
       stat = out$groups
       data$stat = stat[as.character(data$combine__grp__for__statistis_sp), ]$groups

     } else{
		# Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
		# return(Tukey_HSD)
		# Tukey_HSD_table = as.data.frame(Tukey_HSD$combine__grp__for__statistis_sp)
		#print(Tukey_HSD$combine__grp__for__statistis_sp)

		if (length(unique(data$combine__grp__for__statistis_sp)) == 2) {
			Tukey.levels = Tukey_HSD$combine__grp__for__statistis_sp[, 4, drop=F]
		} else {
			Tukey.levels = Tukey_HSD$combine__grp__for__statistis_sp[, 4]
		}

		Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
		Tukey.labels$group = rownames(Tukey.labels)
		Tukey.labels = Tukey.labels[order(Tukey.labels$group), ]
		data$stat = Tukey.labels[as.character(data$combine__grp__for__statistis_sp),]$Letters
    }

    max = max(data[, c(yvariable)])
    min = min(data[, yvariable])
    x = data[, c(xvariable, yvariable, "combine__grp__for__statistis_sp")]
    y = x %>% group_by(combine__grp__for__statistis_sp) %>% summarise(Max =
                                                                        max(!!yvariable_en))
    y = as.data.frame(y)
    # print(y)
    colnames(y) <- c("group", "Max")
    rownames(y) = y$group
    data$y = y[as.character(data$combine__grp__for__statistis_sp),]$Max * 1.04
    # print(data)

    p <- p + geom_text(
      data = data,
      aes(
        x = !!xvariable_en,
        y = y,
        color = !!legend_variable_en,
        label = stat,
        group = !!legend_variable_en
      ),
      position = position_dodge(width =
                                  0.9),
      show.legend = F
    )

    p <- sp_manual_color_ggplot2(p,
                                 data,
                                 legend_variable,
                                 manual_color_vector)
  }

  if (outlier) {
    stats <- boxplot.stats(data[[yvariable]])$stats
    ylim_zoomin <- c(stats[1] / out_scale , stats[5] * out_scale)
    p <- p + coord_cartesian(ylim = ylim_zoomin)
  }

  additional_theme = list()

  if (!sp.is.null(facet_variable)) {
    if (facet_singlecell_style) {
	  if(coordinate_flip){
        strip.position = "top"
	  } else {
	    strip.position = 'left'
	  }
      p <-
        p + facet_wrap(
          ~  .data[[facet_variable]],
          ncol = facet_ncol,
          nrow = facet_nrow,
          scales = facet_scales,
          strip.position = strip.position,
          labeller = as_labeller(unique(data[facet_variable]))
        )
      additional_theme$strip.background = element_blank()
      additional_theme$strip.placement = "outside"
      additional_theme$panel.background = element_rect(fill = NA, colour = 'black')
      y_label = ''
    } else {
      p <-
        sp_ggplot_facet(p, facet_variable, facet_ncol, facet_nrow, facet_scales)
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
    filename = filename,
    additional_theme = additional_theme,
    ...
  )
  p
}
