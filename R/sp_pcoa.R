#' Generating pcoa plot
#'
#' @param data Data file With header line, the first column is the rowname, tab seperated. Each row represents variable (normally genes), each column represents samples.
#' @param grp_file Sample group file with first column as sample names, other columns as sample attributes. Below, color, size, shape variable should be existed in this file. If not supplied, each sample will be treated as one group.
#' And a variable named 'group' can be used to set as color or shape variable.
#' @param shape The variable for point shape.
#' @param shape_order The order for shape variable. Default alphabetical order, accept a string like "'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'".
#' @param color The variable for group information. Necessary, this variable will be used as group information and color variable as also.
#' @param color_order The order for color variable. Default alphabetical order, accept a string like "'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'".
#' @param color_v Manually specified colors. Default system default. Accept string in format like '"green", "red"' (both types of quotes needed).
#' @param title Title of picture.
#' @param size The variable for point size. Optional, such as a number or a variable like count, normally should be number column.
#' @param label Label points. Default no-label (FALSE), accept TRUE.
#' @param label_font_size Label font size. Default system default. Accept numbers less than 5 to shrink fonts.
#' @param ... 
#'
#' @return A ggplot2 object 
#' @export
#'
#' @examples
#' 
#' ## Not run:
#' pcoa_data = "pcoa.data"
#' group_pcoa_data = "group_pcoa.data"
#' sp_pcoa(data=pcoa_data,grp_file = group_pcoa_data, color = "genotype")
#' ## End(Not run)
#' 
sp_pcoa <- function(data,
                    grp_file,
                    shape = 'c_t_c_t0304',
                    shape_order = NULL,
                    color = '',
                    color_order = NULL,
                    color_v = NULL,
                    title = "",
                    size = 'ct___',
                    label = FALSE,
                    label_font_size = NULL,
                    debug = FALSE,
                    file = 'data',
                    ...) {
  if (class(data) == "character") {
    file <- data
    data <- sp_readTable(data, row.names = 1)
  }
  
  
  sampleL = rownames(data)
  
  if (is.null(grp_file)) {
    color = "group"
    grp_file <- data
    grp_file[[color]] = sampleL
  } else {
    if (class(grp_file) == "character") {
      grp_file <- sp_readTable(grp_file, row.names = NULL)
      rownames(grp_file) <- grp_file[, 1]
      grp_file <- grp_file[match(sampleL, rownames(grp_file)),]
    } else {
      grp_file <- grp_file
    }
  }
  
  
  if (length(shape_order) > 1) {
    grp_file[[shape]] <-
      factor(grp_file[[shape]], levels = shape_order, ordered = T)
  }
  
  
  if (length(color_order) > 1) {
    grp_file[[color]] <-
      factor(grp_file[[color]], levels = color_order, ordered = T)
  }
  
  
  if (shape != "c_t_c_t0304") {
    shape_level <- length(unique(grp_file[[shape]]))
    shapes = (1:shape_level) %% 30
  }
  
  # vegan:cmdscale计算矩阵矩阵中主坐标轴坐标，取前3维
  pcoa = cmdscale(data, k = 3, eig = T) # k is dimension, 3 is recommended; eig is eigenvalues
  points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  eig = pcoa$eig
  colnames(points) = c("x", "y", "z")
  points = cbind(points, grp_file[rownames(points), ])
  #points\$label_ct = rownames(points)
  
  
  # plot PCo 1 and 2
  
  color_en = sym(color)
  p = ggplot(points, aes(
    x = x,
    y = y,
    color = !!color_en,
    group = !!color_en
  )) +
    geom_point(alpha = alpha) + coord_fixed() +
    labs(
      x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep =
                  ""),
      y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep =
                  ""),
      title = title
    )  +
    stat_ellipse(level = 0.68, na.rm = TRUE) + theme_classic()
  
  if (length(color_v) >= 2) {
    if (is.numeric(grp_file[[color]])) {
      p <-
        p + scale_colour_gradient(low = color_v[1],
                                  high = color_v[2],
                                  name = color)
    } else {
      p <- p + scale_color_manual(values = color_v)
    }
  }
  
  if (size == "ct___") {
    size = 1
    #p <- p + aes(size=1)
  } else {
    p <- p + aes(size = size)
  }
  #p <- p + aes_(size=size)
  
  if (shape  != "c_t_c_t0304") {
    p <- p + aes(shape = shape)
    if (shape_level > 6) {
      p <- p + scale_shape_manual(values = shapes)
    }
  }
  
  if (label) {
    if (!is.null(label_font_size)) {
      p <- p + geom_text_repel(
        label = paste(rownames(points)),
        show.legend = F,
        size = label_font_size
      )
    } else {
      p <- p + geom_text_repel(label = paste(rownames(points)),
                               show.legend = F)
    }
  }
  
  sampFile = as.data.frame(grp_file[[color]], row.names = row.names(grp_file))
  colnames(sampFile) <- c(color)
  

  # loop for each group pair
  data_table <- as.matrix(data)
  
  compare_data <- as.vector(unique(grp_file[[color]]))
  len_compare_data <- length(compare_data)
  for (i in 1:(len_compare_data - 1)) {
    for (j in (i + 1):len_compare_data) {
      tmp_compare <-
        as.data.frame(cbind(sampA = compare_data[i], sampB = compare_data[j]))
      da_adonis(tmp_compare)
    }
  }
  
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
    
  }
  p
}
