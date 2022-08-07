#' Generating pca plot
#'
#' @param data Data file. With header line, the first column is the rowname, tab separated. Each row represents variable (normally genes), each column represents samples.
#' @param group_data Sample group file with first column as sample names, other columns as sample attributes. Below, color, size, shape variable should be existed in this file.
#' @param title Title of picture. Default empty title
#' @param scale Scale data for PCA analysis. Default, prcomp will centralized data by minus mean value and
#' normalize data by column standard deviation dividing.
#' Often, we would normalize data. Only when we care about the real number changes other than the trends,
#' we do not need scale. When this happens, we expect the ranges of data is small for example log-transformed data.
#' @param minimum_mad Minimum mad to keep. Larger mad, larger variance. Default 0.5.
#' @param top_n Use top n most changed variables for PCA computation. Default 0 meaning use all variables.
#' @param color_variable The variable for point color. Optional, such as color.
#' @param manual_color_vector Manually specified colors. Default system default. Accept string in format like `'"green", "red"'`` (both types of quotes needed).
#' @param color_variable_order The order for color variable. Default alphabetical order, accept a string like "'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'".
#' @param log_transform Log-transform data before principle component analysis. Accept NULL, log2, log10.
#' @param log_add Add a value before log-transform to avoid log(0). Default 0 meaning the minimum non-zero value would be used.
#' @param size_variable The variable for point size. Optional, such as a number or a variable like count, normally should be number column.
#' @param shape_variable The variable for point shape. Optional, such as shape.
#' @param shape_variable_order The order for shape variable. Default alphabetical order, accept a string like "'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'".
#' @param dimensions Dimensions to plot. Default 2.	Accept 3 (both color and shape variable needed and should be same variable).
#' @param alpha Transparency value for points. Optional, such as a number or a variable indicating one data column, normally should be number column.
#' @param label_points Label points (using geom_text_repel). Default no-label (FALSE), accept TRUE.
#' @param label_font_size Label font size. Default system default. Accept numbers less than 5 to shrink fonts.
#' @param coord_fixed_ratio Default 1. Specify 0 to tunr this off. A fixed scale coordinate system forces a specified ratio between the physical representation of data units on the axes. The ratio represents the number of units on the y-axis equivalent to one unit on the x-axis. The default, ratio = 1, ensures that one unit on the x-axis is the same length as one unit on the y-axis. Ratios higher than one make units on the y axis longer than units on the x-axis, and vice versa.
#' @inheritParams sp_manual_color_ggplot2
#' @param ... Other parameters given to \code{\link{sp_ggplot_layout}}.
#'
#' @return pdf and xls files.
#' @export
#'
#' @examples
#' pca_test_data <- matrix(runif(3000,0,100000),ncol=6)
#' colnames(pca_test_data) <- c(paste("wt",1:3,sep = ""),paste("ko",1:3,sep = ""))
#' rownames(pca_test_data) <- c(ID = paste0("ENSG",c(1:500)))
#' pca_data <- as.data.frame(pca_test_data)
#' sp_pca(data = pca_data, group_data = NULL)
#'
#'
#' ## Not run:
#' data = "pca.data"
#' group_data = "pca_group.data"
#' sp_pca(data = data, group_data = group_data, color_variable="Conditions", size = "Diameters", shape_variable = "Batch", label = TRUE)
#'
#' sp_pca(data = data, group_data = group_data, color_variable="Conditions", size = "Diameters", shape_variable = "Batch", label = FALSE,dimensions = 3)
#' ## End(Not run)
#'

sp_pca <- function(data,
                   group_data = NULL,
                   title = NULL,
                   scale = TRUE,
                   color_variable = NULL,
                   manual_color_vector = NULL,
                   log_transform = NULL,
                   facet = NULL,
                   size_variable = NULL,
                   shape_variable = NULL,
                   color_variable_order = NULL,
                   top_n = 0,
                   shape_variable_order = NULL,
                   dimensions = 2,
                   alpha = 1,
                   label_points = FALSE,
                   log_add = 0,
                   label_font_size = NULL,
                   minimum_mad = 0,
                   debug = FALSE,
                   filename = NULL,
                   legend.position = "right",
                   coord_fixed_ratio=1,
                   ...) {


  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if ("character" %in% class(data)) {
    data <- sp_readTable(data, row.names = NULL)
    rownames_data <- make.unique(as.vector(data[, 1]))
    data <- data[, -1, drop = F]
    rownames(data) <- rownames_data
  }

  data <- data[var(data) != 0, ]

  if(minimum_mad + top_n != 0){
    data$mad <- apply(data, 1, mad)
    if(minimum_mad>0){
      data <- data[data$mad > minimum_mad , ]
    }
    data_row_num <- dim(data)[1]
    if (top_n != 0 & top_n < data_row_num) {
      data <-
        data[order(data$mad, decreasing = T),,drop=F]
      data <- data[1:top_n, ]
    }
    data <- data[,-dim(data)[2],drop=F]
  }

  data <- as.data.frame(t(data))

  #print(data[1:3,1:5])

  if (!sp.is.null(log_transform)) {
    # print(y_add)
    # Give the minimum non-zero value to add to avoid log2(0)
    if (log_add == 0) {
      log_add = sp_determine_log_add(data)
    }

    data <- data + log_add
    if (log_transform == "log2") {
      data <- log2(data)
    }else if (log_transform == "log10") {
      data <- log10(data)
    }
  }

  sampleL = rownames(data)


  if (sp.is.null(group_data)) {
    data_t_label <- data
    data_t_label$group = sampleL
    data_t_label$Row.names = sampleL
  } else {
    if (class(group_data) == "character") {
      group_data <- sp_readTable(group_data, row.names = NULL)
      rownames(group_data) <- group_data[, 1]
    }
    #print(colnames(group_data))

    data_t_label <- merge(data, group_data, by = 0, all.x = T)
    rownames(data_t_label) <- data_t_label$Row.names
    data_t_label <-
      data_t_label[match(sampleL, data_t_label$Row.names), ]

    #print(data_t_label[1:4,1:5])
  }

  data_colnames <- colnames(data_t_label)

  #print(data_colnames)

  if (!sp.is.null(shape_variable)) {
    if(! (shape_variable %in% data_colnames )){
      stop(paste(shape_variable,'must be column names of data!'))
    }
    if (!sp.is.null(shape_variable_order)){
      data_t_label = sp_set_factor_order(data_t_label, shape_variable, shape_variable_order)
    }
    shape_level <- length(unique(data_t_label[[shape_variable]]))
    shapes = (1:shape_level) %% 30
  }

  if (!sp.is.null(color_variable)) {
    if(! (color_variable %in% data_colnames )){
      stop(paste(color_variable,'must be column names of data!'))
    }
    if (!sp.is.null(color_variable_order)){
      data_t_label = sp_set_factor_order(data_t_label, color_variable, color_variable_order)
    }
  }


  pca <- prcomp(data, scale = scale)

  rotation = pca$rotation
  x = pca$x
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)

  percentVar2 <- as.data.frame(percentVar)
  rownames(percentVar2) <- colnames(x)

  if(!sp.is.null(filename)){
    sp_writeTable(rotation,
                  file = paste0(filename, ".weights.xls"),
                  keep_rownames = T)
    sp_writeTable(x,
                  file = paste0(filename, ".pcs.xls"),
                  keep_rownames = T)
    sp_writeTable(percentVar2,
                  file = paste0(filename, ".pc_variance.xls"),
                  keep_rownames = T)
  }


  if (dimensions == 2) {
    p <-
      autoplot(pca, data = data_t_label, alpha = alpha, scale=0) +
      ggtitle(title)

    if (!sp.is.null(size_variable)) {
      if(! (size_variable %in% data_colnames )){
        stop(paste(size_variable,'must be column names of data!'))
      }
      size_variable_en = sym(size_variable)
      p <- p + aes(size = !!size_variable_en)
    }
    if (!sp.is.null(color_variable)) {
      color_en = sym(color_variable)
      p <- p + aes(colour = !!color_en)
      p <- sp_manual_color_ggplot2(p,
                                   data,
                                   color_variable,
                                   manual_color_vector)
    }

    if (!sp.is.null(shape_variable)) {
      shape_en = sym(shape_variable)
      p <- p + aes(shape = !!shape_en)
      if (shape_level > 6) {
        p <- p + scale_shape_manual(values = shapes)
      }
    }

    if (label_points) {
      if (!sp.is.null(label_font_size)) {
        p <-
          p + geom_text_repel(aes(label = Row.names),
                              show.legend = F,
                              size = label_font_size)
      } else {
        p <- p + geom_text_repel(aes(label = Row.names), show.legend = F)
      }
    }

    x_label = paste0("PC1 (", round(percentVar[1] * 100,1), "% variance)")
    y_label = paste0("PC2 (", round(percentVar[2] * 100,1), "% variance)")

    if(coord_fixed_ratio>0){
      p <- p + coord_fixed(coord_fixed_ratio)
    }


    p <- sp_ggplot_layout(
      p,
      filename = filename,
      legend.position = legend.position,
      x_label = x_label,
      y_label = y_label,
      title = title,
      ...
    )
    p

  } else {
    library(scatterplot3d)
    if (color_variable != "c_t_c_t0304") {
      group = data_t_label[[color_variable]]
      colorA <- rainbow(length(unique(group)))

      colors <- colorA[as.factor(group)]

      colorl <- colorA[as.factor(unique(group))]
    }

    if (shape_variable != "c_t_c_t0304") {
      group <- data_t_label[[shape_variable]]
      pch_l <- as.numeric(as.factor(unique(group)))
      pch <- pch_l[as.factor(group)]
    }

    pc <- as.data.frame(pca$x)

    saveplot = paste0(filenames, mid, ".pdf")

    if (!sp.is.null(saveplot)) {
      base_plot_save(saveplot, ...)
    }

    # pdf(paste0(filename,mid,"sds.pdf"))
    scatterplot3d(
      x = pc$PC1,
      y = pc$PC2,
      z = pc$PC3,
      pch = pch,
      color = colors,
      xlab = paste0("PC1 (", round(percentVar[1] * 100), "% variance)"),
      ylab = paste0("PC2 (", round(percentVar[2] * 100), "% variance)"),
      zlab = paste0("PC3 (", round(percentVar[3] * 100), "% variance)")
    )

    legend(
      -3,
      8,
      legend = levels(as.factor(color)),
      col = colorl,
      pch = pch_l,
      xpd = T,
      horiz = F,
      ncol = 6
    )

    if (!sp.is.null(saveplot)) {
      dev.off()
    }
  }

}
