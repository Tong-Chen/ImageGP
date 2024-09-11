#' Generating upsetView plot
#'
#' Input file is a matrix:
#'
#' vennFormat 0
#'
#' (First row would be treated as header line. First column is just a normal column (but needed). 0 represents the sample does not contain the genes in row. 1 represents the containing relationship)
#'
#' ID	Samp1	Samp2	Samp3	Samp4	Samp5
#'
#' G1	1	0	1	0	1
#'
#' G2	0	0	1	1	1
#'
#' G3	1	1	1	0	1
#'
#' G4	1	1	1	0	0
#'
#' G5	0	1	0	1	1
#'
#' G6	1	0	1	0	0
#'
#' vennFormat 1 or 2
#'
#' The output contains two barplots, horizontal bar represents the number of genes in each sample, which is the sum of all 1 in sample column. Vertical bar represents the number of sample specific and common genes as indicated by linking vertical lines and points (just as the overlapping regions of venndiagram).
#'
#'
#'
#' @param data Data file. Receive long and wide table forms.
#' @param vennFormat Venn diagram format without header line. Default 0 represents normal data. Accept 1,2.
#' 0: represents wide data listed above.
#' 1: represents venn diagram format without header line.
#' 2: represents venn diagram format with header line.
#' @param sets Specific sets to look at (Include as combinations. Ex: c('Name1', 'Name2')).
#' @param nintersects Number of intersections to plot. If set to NA, all intersections will be plotted.
#' @param order.by How the intersections in the matrix should be ordered by. Options include frequency (entered as 'freq'), degree.
#' @param decreasing How the variables in order.by should be ordered. 'freq' is decreasing (greatest to least) and 'degree' is increasing (least to greatest).
#' @param scale.intersections The scale to be used for the intersection sizes. Options: 'identity', 'log10', 'log2'.
#' @param scale.sets The scale to be used for the set sizes. Options: 'identity', 'log10', 'log2'.
#' @param queries_bar1 Specifies an intersection. Changes the column color.
#' @param queries_bar1_color Input color. Specifies an intersection to use this color.
#' @param pointsize Point size. Default 8.
#' @param maxsets Maximum allowed number of sets. Default 100.
#' @param keep_empty Keep empty intersections. Default FALSE. Accept TRUE to remove empty intersections.
#' @inheritParams base_plot_save
#' @param ... Other parameters given to `base_plot_save`
#'
#' @return A pdf file.
#' @export
#'
#' @examples
#'
#' upsetview_data <- data.frame(elements=c("1","2","2","2","3"), sets=c("A","A","B","C","C"))
#' sp_upsetview(data = upsetview_data, vennFormat=2, saveplot = "upsetView_long.pdf")
#'
#'
#' ## Not run:
#' upsetview_data = "upsetview.data"
#' sp_upsetview(data = upsetview_data, saveplot = "upsetView_wide.pdf")
#' ## End(Not run)
#'
sp_upsetview <- function (data,
                          vennFormat = 0,
                          pointsize = 8,
                          keep_empty = FALSE,
                          sets = NULL,
                          nintersects = NA,
                          order.by = "freq",
                          decreasing = TRUE,
                          scale.intersections = "identity",
                          scale.sets = "identity",
                          queries_bar1 = NULL,
                          queries_bar2 = NULL,
                          queries_bar3 = NULL,
                          queries_bar1_color = NULL,
                          queries_bar2_color = NULL,
                          queries_bar3_color = NULL,
                          saveplot = NULL,
                          debug = FALSE,
                          saveppt = FALSE,
						  maxsets = 100,
                          main_bar_color_vector = "gray23",
                          constantColor =T,
                          ...) {
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if (vennFormat == 0) {
    if ("character" %in% class(data)) {
      data <- sp_readTable(data, row.names = NULL)
    } else if (class(data) != "data.frame") {
      stop("Unknown input format for `data` parameter.")
    }
    data[,-1][data[,-1] != 0] <- 1
    data[,-1][data[,-1] == 0] <- 0

    #long_data <- melt(data)
    #long_data<- long_data[long_data$value == 1,][,1:2]
    #colnames(long_data) <- c("Item", "Grp")
    #sp_writeTable(long_data,file="set_count_input.txt",keep_rownames =  F)
    #system("python set_count.py -f set_count_input.txt -c Item -o jx")
    # source_python("./set_count.py")
    #set_count_data <- sp_readTable("jx.txt", row.names = NULL,header = F)
    #number<- nrow(set_count_data)

  } else {
    header = ifelse(vennFormat == 1, F, T)
    if ("character" %in% class(data)) {
      data <- sp_readTable(data, row.names = NULL, header = header)
    } else if (class(data) != "data.frame") {
      stop("Unknown input format for `data` parameter.")
    }
    data <- unique(data[, 1:2])
    colnames(data) <- c("Item", "Grp")

    #sp_writeTable(data,file="set_count_input.txt",keep_rownames =  F)
    #system("python set_count.py -f set_count_input.txt -c Item -o jx")
    # source_python("./set_count.py")
    #set_count_data <- sp_readTable("jx.txt", row.names = NULL,header = F)
    #number<- nrow(set_count_data)

    data = as.data.frame(reshape2::acast(data, Item ~ Grp, length))
    data = cbind(ID = rownames(data), data)
  }

  data

  nsets = dim(data)[2] - 1

  if(nsets>maxsets){
  	stop(paste0("Your fileformat maybe wrong. It is almost meaningless for interactions among more than ", maxsets, "sets."))
  }

  if (keep_empty) {
    keep_empty = 'on'
  } else {
    keep_empty = NULL
  }

  list1 = list2 = list3 = NULL
  if (!sp.is.null(queries_bar1)) {
    list1 = list(
      query = intersects,
      params = list(queries_bar1),
      color = queries_bar1_color,
      active = T
    )
  }
  if (!sp.is.null(queries_bar2)) {
    list2 = list(
      query = intersects,
      params = list(queries_bar2),
      color = queries_bar2_color,
      active = T
    )
  }
  if (!sp.is.null(queries_bar3)) {
    list3 = list(
      query = intersects,
      params = list(queries_bar3),
      color = queries_bar3_color,
      active = T
    )
  }
  queries_para = list(list1, list2, list3)
  if (!sp.is.null(list1) &&
      !sp.is.null(list1) && !sp.is.null(list1)) {
    queries_para1 = queries_para[!sapply(queries_para, is.null)]
  } else {
    queries_para1 = NULL
  }
  # if (sp.is.null(queries_bar1) && sp.is.null(queries_bar2) && sp.is.null(queries_bar3) ){
  #   queries_para = NULL
  # } else {
  # queries_para =  eval(parse(text = paste(
  #   "list(", list1,")")))
  # }


  if (main_bar_color_vector != "gray23") {
    keep_empty = 'on'
    print(nintersects)

    if (!sp.is.null(keep_empty)){
      all_possible_nintersects = sum(choose(nsets, 1:nsets))

      if (!is.na(nintersects)){
        main_bar_color_vector <-
          c(generate_color_list(main_bar_color_vector, nintersects,
                              alpha = 1, constantColor = T),
            rep('gray23', all_possible_nintersects-nintersects))

      } else {
        main_bar_color_vector <-
          generate_color_list(main_bar_color_vector, all_possible_nintersects,
                              alpha = 1, constantColor = T)

      }
    }else{
      main_bar_color_vector = "gray23"
    }


  }

  a = UpSetR::upset(
    data,
    sets = sets,
    nsets = nsets,
    nintersects = nintersects,
    sets.bar.color = "#56B4E9",
    order.by = order.by,
    decreasing = decreasing,
    scale.intersections = scale.intersections,
    scale.sets = scale.sets,
    empty.intersections = keep_empty,
    queries = queries_para1,
    main.bar.color = main_bar_color_vector
  )


  # list1 = list(list(query=intersects, params=list("Samp1","Samp3"), color="red", active=T))


  if (!sp.is.null(saveplot)) {
    base_plot_save(saveplot, pointsize = pointsize, ...)
    print(a)
    dev.off()
  }
  if (saveppt) {
    p <- UpSetR::upset(
      data,
      sets = sets,
      nsets = nsets,
      nintersects = nintersects,
      sets.bar.color = "#56B4E9",
      order.by = order.by,
      decreasing = decreasing,
      scale.intersections = scale.intersections,
      scale.sets = scale.sets,
      empty.intersections = keep_empty,
      queries = queries_para1,
      main.bar.color = main_bar_color_vector
    )
    eoffice::topptx(p, filename = paste0(saveplot, ".pptx"))
    dev.off()
  }
  a

}

