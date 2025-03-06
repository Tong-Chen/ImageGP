



#' Generating bar plot
#'
#' @param data Data frame or data file (with header line, the first column will
#' not be treated as row names, tab separated).
#' @param melted `TRUE` for dealing with long format matrix, the program will skip melt preprocess. If input is wide format matrix, this parameter should be set to `FALSE`.
#' @param xvariable Name for x-axis variable (one of column names, should be specified
#' when inputting long format matrix).
#' @param color_variable Name for specifying bars colors (one of column names, should be specified
#' when inputting long format matrix).
#' @param color_variable_order Set orders of color variable (this can also used to extract specific rows).
#' @param yvariable Name for value column (one of column names, should be specified
#' when inputting long format matrix).
#' @param xvariable_order Levels for x-axis variable, suitable when x-axis is not used as a number.
#' @param group_variable Specify group info for for computing means and SDs.
#' @param add_point Set to TRUE to add each point (specially used when displaying mean values)
#' @param stat The ways to show the height of bars.
#' The height of bars represent the numerical values in each group by default (normally in `yvariable` column of melted data).
#' One can also give `count` to let the program count the number of
#' items in each group (Normally the `color_variable` column is used to group
#' 'xvariable' column after melt).
#' Or one can give `weight` which will sum values of each group.
#' Default `identity`, accept `count` when categorical data are given.
#' @param bar_mode The ways to place multiple bars for one group.
#' Multiple bars in same place will be stacked together by default.
#' Giving `fill` to get stacked percent bar-plot.
#' Giving `dodge` to arrange multiple bars side-by-side.
#' Default `stack`, accept `dodge`, `fill`.
#' @param facet_variable_order The levels of wrapping to set the order of each group.
#' @inheritParams sp_ggplot_facet
#' @inheritParams sp_transfer_one_column
#' @param error_bar_variable Error-bar column (one of column names). Specify the column containing error bars.
#' @param base_font_size Font-size. Default 11.
#' @param extra_ggplot2_cmd Other legal R codes for ggplot2 will be given here.
#' @param xtics Display xtics. Default TRUE.
#' @param ytics Display ytics. Default FALSE.
#' @param add_text 	Add text to bar. Default FALSE.
#' @inheritParams sp_load_font
#' @inheritParams sp_boxplot
#' @inheritParams sp_ggplot_layout
#' @inheritParams sp_manual_fill_ggplot2
#' @inheritParams dataFilter2
#' @inheritParams sp_read_in_long_wide_matrix
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @param ...
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' bar_test_data <- data.frame(ID = letters[1:4], Gene = letters[c(8,8,9,9,10,10,11,11)], Expr = runif(16))
#' sp_barplot(data = bar_test_data, xvariable = "ID", yvariable = "Expr", color_variable = "Gene")
#'
#' ## Not run:
#' bar_data = "bar.data"
#'
#' sp_barplot(data = bar_data, xvariable = "ID", yvariable = "Expr", color_variable = "Gene")
#' ## End(Not run)
#'
sp_barplot <- function (data,
                        metadata = NULL,
                        color_variable = NULL,
                        yvariable = NULL,
                        xvariable = NULL,
                        melted = TRUE,
                        title = NULL,
                        x_label = NULL,
                        y_label = NULL,
                        top_n = 1,
                        statistical_value_type = sum,
                        keep_filtered_as_others = TRUE,
                        color_variable_order = NULL,
                        xvariable_order = NULL,
                        y_add = 0,
                        group_variable = NULL,
                        add_bar_link = FALSE,
                        add_point = F,
                        yaxis_scale_mode = NULL,
                        facet_variable = NULL,
                        stat = 'identity',
                        bar_mode = 'stack',
                        facet_variable_order = NULL,
                        facet_nrow = NULL,
                        facet_ncol = NULL,
                        error_bar_variable = NULL,
                        base_font_size = 11,
                        legend.position = 'right',
                        xtics = TRUE,
                        xtics_angle = 0,
                        statistics = FALSE,
                        ytics = TRUE,
                        manual_color_vector = "Set2",
                        facet_scales = 'fixed',
                        extra_ggplot2_cmd = NULL,
                        coordinate_flip = FALSE,
                        add_text = FALSE,
                        font_path = NULL,
                        debug = FALSE,
                        filename = NULL,
                        ...) {
  options(scipen = 999)

  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  fontname = sp_load_font(font_path = font_path)


  if (melted) {
    if (sp.is.null(xvariable) || sp.is.null(yvariable)) {
      stop("For melted matrix, <xvariable> and <yvariable> should be supplied.")
    }
  } else {
      #if (sp.is.null(color_variable)) {
        color_variable = "color_variable"
      #}
      xvariable = 'variable'
      yvariable = "value"
  }

  #print(xvariable)

  data <- sp_read_in_long_wide_matrix(data, color_variable, melted,
                                      top_n = top_n,
                                      statistical_value_type = statistical_value_type,
                                      keep_filtered_as_others = keep_filtered_as_others)

  #print(data)

  wide_rownames <- data$wide_rownames
  wide_colnames <- data$wide_colnames
  data <- data$data

  if (!sp.is.null(metadata) && metadata != "") {
    if (class(metadata) == "character") {
      metadata <- sp_readTable(metadata, row.names = NULL)
    } else if (!"data.frame" %in% class(data)) {
      stop("Unknown input format for `metadata` parameter.")
    }

    # matched_column <-
    #   get_matched_columns_based_on_value(data, metadata,
    #                                      only_allow_one_match =
    #                                        T)
    #
    # # return(list(data=data, metadata=metadata, matched_column=matched_column))
    # data <-
    #   merge(data, metadata, by.x = matched_column[1], by.y = matched_column[2], suffixes = c("",".y"))
    # data[[matched_column[2]]] = data[[matched_column[1]]]
    data <- merge_data_with_auto_matched_column(data, metadata)
  }

  # print(data)

  data_colnames <- colnames(data)

  if (!(xvariable %in% data_colnames &&
        yvariable %in% data_colnames)) {
    stop(paste(xvariable, 'or', yvariable, 'must be one of column names of data!'))
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



  #print(data)

  xvariable_en = sym(xvariable)

  yvariable_en = sym(yvariable)

  point_yvariable_en = yvariable_en

  data_point = data

  if(!sp.is.null(group_variable)){
    if (!(group_variable %in% data_colnames)) {
      stop(paste(group_variable,'must be one of column names of data!'))
    }
    # group_variable_en = sym(group_variable)
    group_variable_vector <- unique(c(xvariable, group_variable, facet_variable))
    group_variable_vector <- group_variable_vector[!sapply(group_variable_vector, sp.is.null)]
    if (length(group_variable_vector) == 1 ){
      xvariable = group_variable_vector
      color_variable = group_variable_vector
    } else {
      color_variable = group_variable
    }
    data_sd_mean <- data %>% group_by(across(group_variable_vector)) %>%
      summarise(Standard_deviation=sd(!!yvariable_en), Mean_value=mean(!!yvariable_en)) %>%
      ungroup() %>%
      group_by(!!xvariable_en) %>%
      mutate(Mean_value_cumsum_s_p=rev(cumsum(rev(Mean_value))))

    data <- as.data.frame(data_sd_mean)
    print(data_sd_mean)

    # data_sd_mean = sp_set_factor_order(data_sd_mean, group_variable, group_variable_order)


    # data <- merge(data, data_sd_mean, by=group_variable, all=F)

    yvariable = "Mean_value"
    yvariable_en = sym(yvariable)

    error_bar_variable = "Standard_deviation"
    error_bar_variable_en = sym(error_bar_variable)

    if(sp.is.null(color_variable)){
      color_variable <- group_variable[group_variable!=xvariable][1]
    }

    #bar_mode = "dodge"
    #print(data)
  }




  if (!melted){
    xvariable_order = wide_colnames
    color_variable_order = wide_rownames
  }

  data = sp_set_factor_order(data, xvariable, xvariable_order)

  #print(data)

  if (!sp.is.null(color_variable) && color_variable != xvariable) {
    if (!(color_variable %in% data_colnames)) {
      stop(paste(color_variable,'must be one of column names of data!'))
    }
    data = sp_set_factor_order(data, color_variable, color_variable_order)
  } else {
    color_variable = xvariable
  }

  color_variable_en = sym(color_variable)

  #print(data)

  if (!sp.is.null(facet_variable)) {
    if (!(facet_variable %in% data_colnames)) {
      stop(paste(facet_variable,'must be one of column names of data!'))
    }
    data = sp_set_factor_order(data, facet_variable, facet_variable_order)
  }

  if (bar_mode  == "fill" && add_text) {
    data <-
      data %>% group_by(!!xvariable_en) %>%
      mutate(count = sum(!!yvariable_en)) %>%
      mutate(freq = round(100 * !!yvariable_en / count, 2))
  }

  if(bar_mode == "stack" && (!"Mean_value_cumsum_s_p" %in% colnames(data))){
    # print(data[[xvariable]])
    data <- data %>% group_by(!!xvariable_en) %>%
      mutate(Mean_value_cumsum_s_p=rev(cumsum(rev(!!yvariable_en))))
    # print(data)
  }


  xvariable_en = sym(xvariable)
  color_variable_en = sym(color_variable)
  yvariable_en = sym(yvariable)

  width_dodge = 0.75
  #print(data)

  if (bar_mode  == "dodge") {
    position = position_dodge(width = width_dodge)
    errorbar_base_variable = yvariable
  }else if (bar_mode  == "stack") {
    position = position_stack(vjust = 0.5)
    errorbar_base_variable = "Mean_value_cumsum_s_p"
  }else if (bar_mode  == "fill") {
    position = position_fill(vjust = 0.5)
    errorbar_base_variable = "Mean_value_cumsum_s_p"
  }

  if (stat == "count") {
    p <- ggplot(data, aes(x = !!xvariable_en, group = !!yvariable_en))
  } else {
    p <-
      ggplot(data,
             aes(
               x = !!xvariable_en,
               y = !!yvariable_en,
               group = !!color_variable_en
             ))
  }

  p <-
    p + geom_bar(
      stat = stat ,
      position = bar_mode ,
      aes(fill = !!color_variable_en),
      width = width_dodge
    )
  data_link<- sp_set_factor_order(data_point, xvariable, xvariable_order)
  if (add_bar_link && bar_mode != "dodge") {
    wild_data <- spread(  data = data_link,  key = xvariable, value = yvariable )
    xvariable_order_link <- as.character(unique(data_link[,xvariable]))
    color_variable_order_link <- as.character(unique(data_link[,color_variable]))
    wild_data[[color_variable]] <- factor(wild_data[[color_variable]],
                               levels = color_variable_order_link, ordered = T)
    wild_data <- wild_data[order(wild_data[,color_variable],decreasing=T),]
    wild_data <- wild_data[, c(color_variable,xvariable_order_link)]
    wild_data_col <- colnames(wild_data)
    wild_data_row <- rownames(wild_data)
    if (sp.is.null(color_variable_order)){
    if (bar_mode == "stack") {
      link_dat <- wild_data %>%
        arrange(by = desc(color_variable)) %>%
        mutate_if(is.numeric, cumsum)
    } else {
      wild_data_colorvariable <- wild_data[color_variable]
      wild_data <-
        cbind(wild_data_colorvariable, as.data.frame(apply(wild_data[, -1], 2, function(x)
          x / sum(x))))
      link_dat <- wild_data  %>%
        arrange(by = desc(color_variable)) %>%
        mutate_if(is.numeric, cumsum)
    }
    } else {
      wild_data = sp_set_factor_order(wild_data, color_variable, color_variable_order)
      wild_data <- wild_data[order(wild_data[,color_variable],decreasing=T),]
      if (bar_mode == "stack") {


        link_dat <- wild_data %>%
          arrange(by = desc(color_variable)) %>%
          mutate_if(is.numeric, cumsum)
      } else {

        wild_data_colorvariable <- wild_data[color_variable]
        wild_data <-
          cbind(wild_data_colorvariable, as.data.frame(apply(wild_data[, -1], 2, function(x)
            x / sum(x))))

        link_dat <- wild_data  %>%
          arrange(by = desc(color_variable)) %>%
          mutate_if(is.numeric, cumsum)
      }
    }

    if (ncol(link_dat) < 4){
      link_dat <- data.frame(y=t(matrix(t(link_dat[,-1]), nrow=2)))

      link_dat$x.1 <- 1:(ncol(wild_data) - 2) + width_dodge / 2
      link_dat$x.2 <- 1:(ncol(wild_data) - 2) + (1 - width_dodge / 2)
      p <- p + geom_segment(data=link_dat, aes(x=x.1, xend=x.2, y=y.1, yend=y.2), inherit.aes = F)
    } else {
    link_dat <-
      link_dat[, c(1, 2, rep(3:(ncol(link_dat) - 1), each = 2), ncol(link_dat))]
    link_dat <- data.frame(y = t(matrix(t(link_dat[, -1]), nrow = 2)))
    link_dat$x.1 <- 1:(ncol(wild_data) - 2) + width_dodge / 2
    link_dat$x.2 <- 1:(ncol(wild_data) - 2) + (1 - width_dodge / 2)

    p <- p + geom_segment(data = link_dat,
                          aes(
                            x = x.1,
                            xend = x.2,
                            y = y.1,
                            yend = y.2
                          ),
                          inherit.aes = F)
    }
  }

  if (!sp.is.null(error_bar_variable)) {
    if (!(error_bar_variable %in% c(data_colnames, "Standard_deviation"))) {
      stop(paste(error_bar_variable,'must be column names of data!'))
    }


    if(bar_mode == "fill"){
      bar_mode = "stack"
   }
    error_bar_variable_en = sym(error_bar_variable)
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
          # position = "identity"
          position = position
        )
    }
  }


  if (bar_mode  == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }

  if (add_point){
    p <- p + geom_quasirandom(aes(x = !!xvariable_en,
                                  y = !!point_yvariable_en,
                                  group=!!color_variable_en),
                              data = data,
                              color = "grey",
                              varwidth = T,
                              groupOnX = TRUE,
                              dodge.width = width_dodge,
                              position =  position)
  }


  if(add_text){
    text_size =  base_font_size / 3.2
    geom_text_parameter <- list()



    geom_text_parameter$position = position

    if(!sp.is.null(fontname)){
      geom_text_parameter$famliy = fontname
    }

    geom_text_parameter$size = text_size
    geom_text_parameter$show.legend = F

    if(sp.is.null(error_bar_variable)){
      sp_geom_text <- function(...){
        ggplot2::geom_text(mapping=aes(label = !!yvariable_en), ...)
      }
      p <-
        p + do.call(sp_geom_text, c(geom_text_parameter))
    } else {
      sp_geom_text1 <- function(...){
        geom_text(mapping=aes(
          label = sprintf("%.2f", !!yvariable_en - !!error_bar_variable_en),
          y = !!yvariable_en - !!error_bar_variable_en
        ),
        vjust = 1.5, ...)
      }
      sp_geom_text2 <- function(...){
        geom_text(mapping=aes(
          label = sprintf("%.2f", !!yvariable_en + !!error_bar_variable_en),
          y = !!yvariable_en + !!error_bar_variable_en
        ),
        vjust = .5, ...)
      }
      p <-
        p + do.call(sp_geom_text1, c(geom_text_parameter)) +
        do.call(sp_geom_text2, c(geom_text_parameter))
    }
  }

  if (statistics) {
    # 代码修改自 amplicon包 microbiota/amplicon
    # https://github.com/microbiota/amplicon/blob/master/R/alpha_boxplot.R

    group_variable_vector <- unique(c(xvariable, color_variable, facet_variable))
    group_variable_vector <- group_variable_vector[!sapply(group_variable_vector, sp.is.null)]
    #data2 <- data[,group_variable_vector]
    data$combine__grp__for__statistis_sp <- do.call(paste0, data[group_variable_vector])

    formula = as.formula(paste(yvariable, "~", "combine__grp__for__statistis_sp"))
    model = aov(formula, data = data)
    if (length(unique(data$combine__grp__for__statistis_sp)) == 2) {
      library(agricolae)
      out = LSD.test(model, "combine__grp__for__statistis_sp", p.adj = "none")
      # print(out)
      LSD.test_table = as.data.frame(out$statistics)
      stats = out$groups
      data$stats = stats[as.character(data$combine__grp__for__statistis_sp), ]$groups

      suppressWarnings(write.table(
        LSD.test_table,
        file = "barplot_LSD.test.txt",
        sep = "\t",
        quote = F,
        row.names = F
      ))
    } else{
      Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
      # return(Tukey_HSD)
      Tukey_HSD_table = as.data.frame(Tukey_HSD$combine__grp__for__statistis_sp)
      Tukey.levels = Tukey_HSD$combine__grp__for__statistis_sp[, 4]
      Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
      Tukey.labels$group = rownames(Tukey.labels)
      Tukey.labels = Tukey.labels[order(Tukey.labels$group),]
      data$stats = Tukey.labels[as.character(data$combine__grp__for__statistis_sp), ]$Letters
      # print(data)
      suppressWarnings(write.table(
        Tukey_HSD_table,
        file = "barplot_TukeyHSD.txt",
        sep = "\t",
        quote = F,
        row.names = F
      ))
    }

    max = max(data[, c(yvariable)], na.rm=T)
    min = min(data[, yvariable], na.rm=T)
    x = data[, c(xvariable, yvariable, "combine__grp__for__statistis_sp")]
    y = x %>% group_by(combine__grp__for__statistis_sp) %>% summarise(Max =
                                                                        max(!!yvariable_en))
    y = as.data.frame(y)
    # print(y)
    colnames(y) <- c("group", "Max")
    rownames(y) = y$group
    data$y = y[as.character(data$combine__grp__for__statistis_sp),]$Max * 1.04
    # print(data)
    p + geom_text(
      data = data,
      aes(
        x = !!xvariable_en,
        y = y,
        color = !!color_variable_en,
        label = stats,
        group = !!color_variable_en
      ),
      position = position_dodge(width =
                                  0.9),
      show.legend = F
    )

    p <- sp_manual_color_ggplot2(p,
                                 data,
                                 color_variable,
                                 manual_color_vector)
  }

  if (!sp.is.null(facet_variable)) {
    p <-
      sp_ggplot_facet(p, facet_variable, facet_ncol, facet_nrow, facet_scales)
  }

  if (!sp.is.null(yaxis_scale_mode) &&
      (yaxis_scale_mode  != "log2") &&
      (yaxis_scale_mode  != "log10")) {
    p <- p +  eval(parse(text = yaxis_scale_mode))
  }


  p <-
    sp_manual_fill_ggplot2(p, data, color_variable, manual_color_vector)



  additional_theme <- list()

  if (!xtics) {
    additional_theme$axis.text.x = element_blank()
  }
  if (!ytics) {
    additional_theme$axis.text.y = element_blank()
  }

  additional_theme$axis.ticks.x = element_blank()
  additional_theme$legend.key  = element_blank()

  if(xvariable == "variable"){
    x_label = " "
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
      coordinate_flip = coordinate_flip,
      additional_theme = additional_theme,
      fontname = fontname,
      base_font_size = base_font_size,
      ...
    )

  p


}
