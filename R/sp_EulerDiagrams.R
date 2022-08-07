#' Generate Euler diagrams (proposional)
#'
#' When `format` is `items`, `data` examples are
#' ```
#' Items  Group
#' g1	Set1
#' g2	Set1
#' a1	Set3
#' a3	Set1
#' b4	Set1
#' g1	Set2
#' h1	Set4
#' ```
#'
#' When `format` is `counts`, `data` examples are
#'
#' ```
#' Intersection	Count
#' Set1&Set2	2
#' Set1&Set3&Set4&Set5	1
#' Set3&Set5	2
#' Set1	1
#' Set2&Set4	1
#' Set2&Set3&Set4	1
#' Set2&Set5	1
#' Set3	1
#' Set5	1
#' Set4&Set5	1
#' Set4	2
#' ```
#'
#' @param data One filename or dataframe containing data in specified formats with header line.
#' @param format `items` or `counts` with format specified above.
#' @param intersection_variable Only used when `format=counts` to specify which column
#' contains different types of interactions. For example data, first column name `Intersection`
#' should be given here. Color should be specified for the appearance order of each set.
#' @param count_variable Only used when `format=counts` to specify which column
#' contains computed counts for different types of interactions. For example data,
#' first column name `Count` should be given here.
#' @param item_variable Only used when `format=items` to specify which column
#' contains all items like genes or OTUs or species. For example data, first column
#' name `Items` should be given here.
#' @param set_variable Only used when `format=items` to specify which column
#' contains group information of items. For example data, second column
#' name `Group` should be given here. Color should be specified for the alphabetial
#' order of each set.
#' @param type Show `percent` or `counts` in the plot. Default `counts`.
#' @param shape Use `circle` or `ellipse` in the plot. Default `circle`.
#' @inheritParams generate_color_list
#' @inheritParams sp_ggplot_layout
#' @param font_quantities Font size for numbers in Euler plot. Default `1`.
#' @param lty Line type of circle or ellipse edges from `1` to `6` represents `solid`,
#' `dashed`, `dotted`, `dotdash`, `longdash` and `twodash` separately. Default `1`.
#' @param labels_font Font size for labels in Euler plot. Default `1`.
#' @inheritParams base_plot_save
#' @param ... Other parameters given to base_plot_save
#'
#' @return a grid object
#' @export
#'
#' @examples
#'
#'
#'
#'
sp_EulerDiagrams <- function (data,
                              format = "items", # items or counts
                              intersection_variable = NULL,
                              count_variable = NULL,
                              item_variable = NULL,
                              set_variable = NULL,
                              type = c("counts"), #percent
                              shape = "circle", #ellipse
                              manual_color_vector = NULL,
                              alpha = 1,
                              legend.position = "right",
                              font_quantities = 1,
                              lty = 1,
                              labels_font = 1,
                              saveplot = NULL,
                              saveppt = FALSE,
                              ...) {

  # library(eulerr)
  # 最基本的输入格式判断不能落下
  if ("character" %in% class(data)) {
    data <- sp_readTable(data, row.names = NULL, header=T)
  } else if (class(data) != "data.frame") {
    stop("Unknown input format for `data` parameter.")
  }

  if(format == "items"){
    if(sp.is.null(item_variable) || sp.is.null(set_variable)){
      stop(paste(item_variable,"and",set_variable,"are required for <items> format data."))
    }
    data <- unique(data[c(item_variable, set_variable)])
    data = split(data[,1], data[,2])
    setname = names(data)
  } else if (format == "counts"){
    if(sp.is.null(intersection_variable) || sp.is.null(count_variable)){
      stop(paste(intersection_variable,"and",count_variable,"are required for <items> format data."))
    }
    #data <- sort(data)
    setname = unique(unlist(strsplit(as.vector(data[[intersection_variable]]),'&')))
    intersection_count <- data[[count_variable]]
    names(intersection_count) = data[[intersection_variable]]
    data = intersection_count
    # print(data)
  } else {
    stop("Only <items> or <counts> format allowed!")
  }

  setnumber <- length(setname)
  if(setnumber>6){
	stop("Datamatrix with more than 6 sets is not allowed!")
  }

  Euler <- euler(data, shape = shape)

  if (!sp.is.null(manual_color_vector)) {
    fill_color <-  generate_color_list(manual_color_vector, setnumber,
                                       alpha = alpha)
    # fill color would be the appearance order of sets in file
    #names(fill_color) <- sort(setname)
    #fill_color <- as.list(fill_color)
  } else {
    fill_color = list()
  }

  print(setname)
  # print(sort(setname))
  # print(fill_color)

  # plot.euler

  a <- plot(
    Euler,
    quantities = list(type = type, font = font_quantities),
    fill = fill_color,
    legend = list(side = legend.position),
    edges = list(lty = lty),
    labels = list(font = labels_font)
  )


  if (!is.null(saveplot)) {
    base_plot_save(saveplot, ...)
    print(a)
    dev.off()
  }

  if (saveppt){
    p<-plot(
      Euler,
      quantities = list(type = type, font = font_quantities),
      fill = fill_color,
      legend = list(side = legend.position),
      edges = list(lty = lty),
      labels = list(font = labels_font)
    )
    eoffice::topptx(p,filename = paste0(saveplot,".pptx"))
    # dev.off()
  }
  a
}
