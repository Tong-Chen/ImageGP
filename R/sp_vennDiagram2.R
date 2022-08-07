#' Generating venDiagram plot
#'
#' @param data Data file (the first column is the name of genes or other things one want to compare, the second column is set name of each items, tab separated).
#' @inheritParams sp_readTable
#' @param supplyNumbers If you have the size of each set and the number overlaps between or among each set,  please give TRUE here.
#' @param title Title for the picture.
#' @param item_variable Specify the column containing all items (one of column names of data).
#' @param set_variable Specify the column containing set names (one of column names of data).
#' @param set_variable_order Specify which sets to be plotted and their order (one of column names of data).
#' (normally whole or a subset of unique values of column `set_variable` with order specified)
#' @param numVector List of numbers for venn plot(used when `supplyNumbers` is true).
#' For two-set venn, the format is "100, 110, 50" represents
#' (length_a, length_b, a_b_overlap).
#' For three-set venn, the format is "100, 110, 90, 50, 40, 40, 20"
#' represents (length_a, length_b, length_c,
#'            a_b_overlap,  b_c_overlap, a_c_overlap, a_b_c_overlap).
#' For four-set venn, the format is "100, 110, 90, 50, 40, 40, 20"
#' represents (length_a, length_b, length_c,
#'            a_b_overlap, a_c_overlap, a_d_overlap, b_c_overlap,
#'            b_d_overlap, c_d_overlap, abc_overlap, abd_overlap,
#'            acd_overlap, bcd_overlap, abcd_overlap).
#' @param labelVector List of label for venn plot(used when `supplyNumbers` is true).
#' Format: c('a', 'b')" for two-set and c('a', 'b', 'c') for three-set.
#' @param manual_color_vector Color for each area. Ussally the number of colors should
#' be equal to the number of labels (however the program will make them equal).
#' If you manually set colors for 4-way
#' venn diagram, the first color will be given to the
#' leftmost set, the second will be given to the rightmost
#' set, the third will be given to second leftmost and the forth
#' will be given to the second rightmost.
#' Colors like c('red', 'blue', '#6181BD') (number of colors not matter) or
#' a RColorBrewer color set like  "BrBG"     "PiYG"     "PRGn"     "PuOr"
#' "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"  "Spectral" "Accent"
#' "Dark2"    "Paired"   "Pastel1"  "Pastel2"  "Set1"
#' "Set2"    "Set3"     "Blues"    "BuGn"     "BuPu"
#' "GnBu"     "Greens"   "Greys"    "Oranges" "OrRd"     "PuBu"
#' "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"
#' "YlGn"    "YlGnBu"   "YlOrBr"   "YlOrRd"
#' (check http://www.sthda.com/english/wiki/colors-in-r for more).
#' @inheritParams sp_manual_color_ggplot2
#' @param label_size Siez of category names. Default system default.
#' @param margin Number giving the amount of whitespace around the diagram in grid units. Default system default
#' @inheritParams base_plot_save
#' @param ...
#'
#' @return A grid object.
#' @export
#'
#' @examples
#'
#'
#' vennDiagram_test_data <- data.frame(elements=c("1","2","2","2","3"), sets=c("A","A","B","C","C"))
#' sp_vennDiagram(data = vennDiagram_test_data, label1 = "A",label2 = "B", label3 = "C")
#'
#'
#'
#' sp_vennDiagram( supplyNumbers = TRUE,  numVector=c (120, 110, 50), labelVector=c('a','b'))
#'
#'
#'
#' ## Not run:
#' vennDiagram_data = "vennDiagram.data"
#' sp_vennDiagram2(data = vennDiagram_data, item_variable = "Gene", set_variable = "Sample",
#' set_variable_order = c("Set1","Set2","Set3", "Set4"))
#' ## End(Not run)
#'
sp_vennDiagram2 <- function (data,
                             supplyNumbers = FALSE,
                             header = TRUE,
                             title = '',
                             item_variable = NULL,
                             set_variable = NULL,
                             set_variable_order = NULL,
                             color_for_circumference = "transparent",
                             numVector = c(),
                             labelVector = c(),
                             manual_color_vector = c("dodgerblue",
                                                     "goldenrod1",
                                                     "darkorange1",
                                                     "seagreen3",
                                                     "orchid3"),
                             alpha = 0.5,
                             label_size = NULL,
                             margin = NULL,
                             filename = NULL,
                             debug = FALSE,
                             saveppt = FALSE,
                             ...) {
  options(stringsAsFactors = FALSE)

  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if (!sp.is.null(filename)) {
    base_plot_save(filename, ...)
  }

  if (!supplyNumbers)  {
    if ("character" %in% class(data)) {
      data <- sp_readTable(data, row.names = NULL, header = header)
    } else if (!"data.frame" %in% class(data)) {
      stop("Unknown input format for `data` parameter.")
    }

    data_colnames <- colnames(data)
    if (sp.is.null(item_variable) || sp.is.null(set_variable)) {
      stop("<item_variable> and <set_variable> should be supplied.")
    }

    if (!(item_variable %in% data_colnames &&
          set_variable %in% data_colnames)) {
      stop(paste(
        set_variable,
        'or',
        item_variable,
        'must be one of column names of data!'
      ))
    }

    data <-
      data[!is.na(data[[item_variable]]), c(item_variable, set_variable)]
    data[[set_variable]] <- make.names(data[[set_variable]])
    data <- unique(data)



    labelRecord = set_variable_order
    labelRecord = labelRecord[labelRecord != "NULL"]

    if (length(labelRecord) == 0) {
      labelRecord = unique(data[[set_variable]])
      if (length(labelRecord) > 5) {
        labelRecord = labelRecord[1:5]
      }
    }

    num = length(labelRecord)

    generateVennData <- function(label, data) {
      label = make.names(label)
      label_content <-
        data[data[[set_variable]] == label, item_variable]
    }

    labelContent = lapply(labelRecord, generateVennData, data)

    names(labelContent) <- labelRecord

    # 移除log
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

    color_v <- generate_color_list(manual_color_vector, num)
    #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
    if (color_for_circumference != "transparent") {
      color_for_circumference <-
        generate_color_list(color_for_circumference, num)
    }

    main.pos = c(0.5, 1.05)
    cat.default.pos = "outer"
    cat.pos = 0
    cat.dist = 0.03
    if (sp.is.null(label_size)) {
      label_size = 1
    }
    if (sp.is.null(margin)) {
      margin = 0.01
    }
    if (num == 5) {
      cat.pos = c(0, 225, 270, 135, 90)
      cat.dist = 0.12
      if (sp.is.null(label_size)) {
        label_size = 0.6
      }
      main.pos = c(0.5, 1)
      if (sp.is.null(margin)) {
        margin = 0.05
      }
    } else if (num == 4) {
      cat.pos = c(-15, 15, 0, 0)
      cat.dist = c(0.22, 0.22, 0.11, 0.11)
    } else if (num == 3) {
      cat.pos = c(0, 0, 180)
    }
    p <- venn.diagram(
      x = labelContent,
      filename = NULL,
      col = color_for_circumference,
      lwd = 1,
      fill = color_v,
      alpha = alpha,
      main = title,
      label.col = c("black"),
      cex = 1,
      cat.col = color_v,
      cat.cex = label_size,
      margin = margin,
      main.pos = main.pos,
      cat.default.pos = cat.default.pos,
      cat.pos = cat.pos,
      cat.dist = cat.dist
    )
  } else {
    #---venn plot for given numbers---------
    num <- length(labelVector)
    color_v <- generate_color_list(manual_color_vector, num)
    #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
    if (color_for_circumference != "transparent") {
      color_for_circumference <-
        generate_color_list(color_for_circumference, num)
    }

    if (num == 2) {
      p <- draw.pairwise.venn(
        area1 = numVector[1],
        area2 = numVector[2],
        cross.area = numVector[3],
        category = labelVector,
        lwd = rep(1, 1),
        lty = rep(2, 2),
        col = color_for_circumference,
        fill = color_v,
        cat.col = color_v,
        ind = F
      )
    } else if (num == 3) {
      p <- draw.triple.venn(
        area1 = numVector[1],
        area2 = numVector[2],
        area3 = numVector[3],
        n12 = numVector[4],
        n23 = numVector[5],
        n13 = numVector[6],
        n123 = numVector[7],
        category = labelVector,
        col = color_for_circumference,
        fill = color_v,
        cat.col = color_v,
        reverse = FALSE,
        ind = F
      )
    } else if (num == 4) {
      p <- draw.quad.venn(
        area1 = numVector[1],
        area2 = numVector[2],
        area3 = numVector[3],
        area4 = numVector[4],
        n12 = numVector[5],
        n13 = numVector[6],
        n14 = numVector[7],
        n23 = numVector[8],
        n24 = numVector[9],
        n34 = numVector[10],
        n123 = numVector[11],
        n124 = numVector[12],
        n134 = numVector[13],
        n234 = numVector[14],
        n1234 = numVector[15],
        category = labelVector,
        col = color_for_circumference,
        fill = color_v,
        cat.col = color_v,
        reverse = FALSE,
        ind = F
      )
    } else if (num == 5) {
      stop("Unsupported 5 set venn diagram for input numbers.")
    }
  }

  if (!sp.is.null(filename)) {
    grid::grid.newpage()
    grid::grid.draw(p)
    dev.off()
  }
  grid::grid.newpage()
  grid::grid.draw(p)

  if (!sp.is.null(filename)) {
    if (saveppt) {
      if (!supplyNumbers)  {
        venn.diagram(
          x = labelContent,
          filename = NULL,
          col = color_for_circumference,
          lwd = 1,
          fill = color_v,
          alpha = alpha,
          main = title,
          label.col = c("black"),
          cex = 1,
          cat.col = color_v,
          cat.cex = label_size,
          margin = margin,
          main.pos = main.pos,
          cat.default.pos = cat.default.pos,
          cat.pos = cat.pos,
          cat.dist = cat.dist
        )
        eoffice::topptx(filename = paste0(filename, ".pptx"))
      }
    }
  }
}
