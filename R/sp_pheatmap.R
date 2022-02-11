
#' Pheatmap function only for inner usages
#'
#' @param n Nothing
#' @param gaps Nothing
#' @param m Nothing
#'
#' @return A list
#' @export
#'
#' @examples
#'
#' #Ignore
#'
find_coordinates <- function(n, gaps, m = 1:n) {
  if (length(gaps) == 0) {
    return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc")))
  }

  if (max(gaps) > n) {
    stop("Gaps do not match matrix size")
  }

  size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))

  gaps2 = apply(sapply(gaps, function(gap, x) {
    x > gap
  }, m), 1, sum)
  coord = m * size + (gaps2 * unit("4", "bigpts"))

  return(list(coord = coord, size = size))
}

#' Pheatmap function only for inner usages
#'
#' @param coln Nothing
#' @param gaps Nothing
#' @param xtics_angle Nothing
#' @param ... Nothing
#'
#' @return A grob
#' @export
#'
#' @examples
#'
#' #Ignore
#'
draw_colnames_custom <-
  function (coln, gaps, xtics_angle = 0, ...) {
    coord = find_coordinates(length(coln),  gaps)
    x = coord$coord - 0.5 * coord$size

    vjust <- 0.5
    hjust <- 0.5

    if (xtics_angle == 90) {
      hjust <- 1
      vjust <- 0
    } else if (xtics_angle >= 180) {
      hjust <- 0.5
      vjust <- -0.5
    } else if (xtics_angle >= 200) {
      hjust <- 0.5
      vjust <- -1
    } else if (xtics_angle == 250) {
      hjust <- 0.5
      vjust <- 0
    } else if (xtics_angle == 0) {
      vjust <- 1
      hjust <- 0.5
    } else {
      vjust <- 1
      hjust <- 1
    }
    #else if (xtics_angle == 90) {
    #  hjust <- 1
    #  vjust <- 0.5
    #} else if (xtics_angle == 0) {
    #  vjust <- 1
    #  hjust <- 0.5
    #}

    res = grid::textGrob(
      coln,
      x = x,
      y = unit(1, "npc") - unit(3, "bigpts"),
      vjust = vjust,
      hjust = hjust,
      rot = xtics_angle,
      gp = gpar(...)
    )
    return(res)
  }



#' Generating pheatmap plot
#'
#' @param data Data file or dataframe (with header line, the first column is the rowname, tab seperated.
#' Colnames normally should be unique unless you know what you are doing.)
#' @param filename Filename for output files.
#' @param renameDuplicateRowNames Specify the way to deal with duplicate row names.
#' Default FALSE: representing duplicated row names are not allowed.
#' Accept  TRUE: representing make duplicated row names unique by adding <.1>, <.2>
#' for the second, third appearances.
#' @param logv First get log-value, then do other analysis.
#' Accept an R function log2 or log10. Default FALSE.
#' @param log_add A value to add before log-transfer in-case log zero.
#' Default 0 the program will automatically choose value to add.
#' @param scale Scale the data or not for clustering and visualization.
#' Default 'none' means no scale, accept 'row', 'column' to scale by row or column.
#' @param annotation_row A file or datafrmae to specify row-annotation with first column
#' same as first column of `data`. Default NULL.
#' @param annotation_col A file or datafrmae to specify col-annotation with first column
#' sanme as first row of `data`. Default NULL.
#' @param cluster_rows Hieratical cluster for rows. Default FALSE, accept TRUE.
#' When there are less than 3 rows or more than 5000 rows, this parameter
#' would always be set to FALSE.
#' @param cluster_cols Hieratical cluster for columns. Default FALSE, accept TRUE.
#' When there are less than 3 columns or more than 5000 columns, this parameter
#' would always be set to FALSE.
#' @param anno_cutree_cols Add column tree-cut results as column annotation.
#' @param anno_cutree_rows Add row tree-cut results as row annotation.
#' @param cluster_cols_variable Reorder branch order of clustered columns by given variable. (Test only)
#' @param cluster_rows_variable Reorder branch order of clustered rows by given variable. (Test only)
#' @param remove_cluster_cols_variable_in_annocol Do not show `cluster_cols_variable` in column annotation.
#' @param remove_cluster_rows_variable_in_annorow Do not show `cluster_rows_variable` in row annotation.
#' @param clustering_method Clustering method, Default "complete".
#' Accept "ward.D", "ward.D2","single", "average" (=UPGMA),
#' "mcquitty" (=WPGMA), "median" (=WPGMC) or "centroid" (=UPGMC)
#' @param clustering_distance_rows Clustering distance method for rows.
#' Default 'pearson', accept 'spearman','euclidean', "manhattan", "maximum",
#' "canberra", "binary", "minkowski", "bray", "kulczynski", "jaccard", "gower", "altGower",
#'  "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis". (Some need vegan package)
#' @param clustering_distance_cols Clustering distance method for cols.
#' Default 'pearson', accept 'spearman','euclidean', "manhattan", "maximum",
#' "canberra", "binary", "minkowski", "bray", "kulczynski", "jaccard", "gower", "altGower",
#'  "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis". (Some need vegan package)
#' @param breaks A sequence of numbers that covers the range of values in mat and
#' is one element longer than color vector. Used for mapping values to colors.
#' Useful, if needed to map certain values to certain colors, to certain values.
#' If value is NA then the breaks are calculated automatically. if value is `quantile`, then
#' the breaks would be computed to generate each quantile.
#' @param breaks_mid Mid value for generating breaks when `quantile` is assigned to break.
#' @param breaks_digits Number of digits kept for breaks. Default 2.
#' @param maximum The maximum value one want to keep, any number larger than given value
#' would be taken as this given maximum value. Default Inf, Optional.
#' @param minimum The smallest value one want to keep, any number smaller will be
#' taken as this given minimum value. Default -Inf, Optional.
#' @param correlation_plot First compute the correlation matrix of given `data`, then
#' heatmap correlation data instead of raw data. Default "None", accept "row" or "col" for
#' row correlation or column correlation.
#' @param xtics_angle Rotation angle for x-axis value. Default 0.
#' @inheritParams sp_boxplot
#' @param fontsize Font size. Default 14.
#' @param manual_annotation_colors_sidebar Annotation color. One can only specify color for each column of
#' row-annotatation or col-annotation. For example,
#' 'class' (two values: C1, C2) and group' (two values:G1, G2) are two row-annotations,
#' 'type' (three values, T1, T2, T3) and 'size' (four values, 1, 2, 3, 4)
#' are two col-annoations.
#' Colors can be specified in a string as `'class=c(C1="blue", C2="yellow"), size=c("white", "green"), type=c(T1="pink", T2="black", T3="cyan")'`
#' or a list as `list(class=c(C1="blue", C2="yellow"),size=c("white", "green"))`.
#' In R, one can use colors() function to get names of all available colors.
#' @param kclu Aggregate the rows using kmeans clustering.
#' This is advisable if number of rows is so big that R cannot
#' handle their hierarchical clustering anymore, roughly more than 1000.
#' Instead of showing all the rows separately one can cluster the
#' rows in advance and show only the cluster centers. The number of clusters can be tuned here.
#' Default 'NA' which means no cluster, other positive interger is accepted for executing
#' kmeans cluster, also the parameter represents the number of expected clusters
#' @param ytics Display ytics.
#' @param xtics Display xtics.
#' @param title Title of picture. Default empty title
#' @param width Picture width
#' @param height Picture height
#' @param saveppt Whether to output PPT format. Default false, doesn't output. Accept TRUE, will output ppt file.
#' @inheritParams pheatmap::pheatmap
#' @param ... Other parameters given to \link[pheatmap]{pheatmap}.
#'
#' @return Generate a PDF and TXT file.
#' @export
#'
#' @examples
#' a = c(12,14,17,11,16)
#' b = c(4,20,15,11,9)
#' c = c(5,7,19,8,18)
#' d = c(15,13,11,17,16)
#' e = c(12,19,16,7,9)
#' pheatmap_data = as.data.frame(cbind(a,b,c,d,e))
#' sp_pheatmap(data = pheatmap_data)
#'
#' ## Not run:
#' pheatmap_data = "pheatmap.data"
#' sp_pheatmap(data = pheatmap_data)
#' ## End(Not run)
#'
#'
sp_pheatmap <- function(data,
                        filename = NA,
                        renameDuplicateRowNames = F,
                        logv = NULL,
                        log_add = 0,
                        scale = 'none',
                        annotation_row = NULL,
                        annotation_col = NULL,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        cluster_cols_variable = NULL,
                        cluster_rows_variable = NULL,
                        remove_cluster_cols_variable_in_annocol = FALSE,
                        remove_cluster_rows_variable_in_annorow = FALSE,
                        clustering_method = 'complete',
                        clustering_distance_rows = 'pearson',
                        clustering_distance_cols = 'pearson',
                        breaks = NA,
                        breaks_mid = NULL,
                        breaks_digits = 2,
                        correlation_plot = "None",
                        maximum = Inf,
                        minimum = -Inf,
                        xtics_angle = 0,
                        manual_color_vector = NULL,
                        fontsize = 14,
                        manual_annotation_colors_sidebar = NULL,
                        cutree_cols = NA,
                        cutree_rows = NA,
                        anno_cutree_cols = F,
                        anno_cutree_rows = F,
                        kclu = NA,
                        ytics = TRUE,
                        xtics = TRUE,
                        width = 0,
                        height = 0,
                        title = '',
                        debug = FALSE,
                        saveppt = FALSE,
                        ...) {
  #filename = 'anything'

  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }


  # Overwrite default draw_colnames with your own version
  assignInNamespace(x = "draw_colnames",
                    value = "draw_colnames_custom",
                    ns = asNamespace("pheatmap"))

  if (class(data) == "character") {
    # if (sp.is.null(outputprefix)) {
    #   outputprefix = data
    #   filename = NA
    # }
    data <-
      sp_readTable(data,
                   row.names = 1,
                   renameDuplicateRowNames = renameDuplicateRowNames)
  } else if(class(data) != "data.frame"){
    stop("Unknown input format for `data` parameter.")
  }

  #print(data)
  # if (sp.is.null(outputprefix)) {
  #   outputprefix = "sp_heatmap"
  #   filename = NA
  # }

  # if (!is.na(filename)) {
  #   filename = paste0(outputprefix, '.pdf')
  # }

  # check numerical

  numeric_check = sapply(data, is.numeric)

  non_numeric_col = names(numeric_check[numeric_check == FALSE])

  if (length(non_numeric_col) > 0) {
    stop(paste(non_numeric_col, "contains non-numeric values."))
  }

  if (!sp.is.null(logv)) {
    if (log_add == 0) {
      log_add = sp_determine_log_add(data)
    }
    # Transfer string to R code
    data <- eval(parse(text=logv))(data + log_add)
  }

  if (!sp.is.null(manual_color_vector)) {
    manual_color_vector <- generate_color_list(manual_color_vector, 100)
  } else {
    manual_color_vector <-
      colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  }

  color_length = length(manual_color_vector)

  legend_breaks = NA
  legend_labels = NA
  # Generate quantile breaks
  if (length(breaks) > 1 || !is.na(breaks)) {
    if (length(breaks) == 1 && breaks == "quantile") {
      summary_v <- summary(c(t(data)))
      summary_v[1] <- summary_v[1]
      summary_v[6] <- summary_v[6]
      if (sp.is.null(breaks_mid)) {
        breaks <-
          unique(c(
            seq(summary_v[1], summary_v[2], length = color_length / 4),
            seq(summary_v[2], summary_v[3], length = color_length / 4),
            seq(summary_v[3], summary_v[5], length = color_length / 4),
            seq(summary_v[5], summary_v[6], length = color_length / 4 -
                  1)
          ))
        legend_breaks <- summary_v
      } else {
        breaks_mid <- as.numeric(breaks_mid)
        breaks <- unique(c(
          seq(summary_v[1], breaks_mid,
              length = color_length / 2),
          seq(breaks_mid, summary_v[6], length = color_length / 2 - 1)
        ))
        legend_breaks <- c(summary_v[1], breaks_mid, summary_v[6])
      }
    } else {
      legend_breaks <- breaks
      length_breaks <- length(breaks)
      if (length_breaks < color_length) {
        # break_cnt <- color_length/length_breaks
        manual_color_vector <-
          generate_color_list(c(manual_color_vector[1], manual_color_vector[100]),
                              length_breaks + 1)
      }
    }

    if (breaks_digits) {
      legend_breaks <-
        as.numeric(prettyNum(legend_breaks, digits = breaks_digits))
    }
    legend_labels <- legend_breaks

    # print(breaks)
  }

  if (!sp.is.null(annotation_row)) {
    if (class(annotation_row) == "character") {
      annotation_row <- sp_readTable(annotation_row, row.names = 1)
      annotation_row <-
        annotation_row[match(rownames(data), rownames(annotation_row)), , drop =
                         F]
    }
    if(!sp.is.null(cluster_rows_variable)){
      if (!cluster_rows_variable %in% colnames(annotation_row) ) {
        stop(paste(
          cluster_rows_variable,
          'must be one of column names of row annotation matrix!'
        ))
      }
    }
  } else {
    annotation_row <- NA
  }

  if (!sp.is.null(annotation_col)) {
    if (class(annotation_col) == "character") {
      annotation_col <- sp_readTable(annotation_col, row.names = 1)
      annotation_col <-
        annotation_col[match(colnames(data), rownames(annotation_col)), , drop =
                         F]
    }

    if(!sp.is.null(cluster_cols_variable)){
      if (!cluster_cols_variable %in% colnames(annotation_col) ) {
        stop(paste(
          cluster_cols_variable,
          'must be one of column names of column annotation matrix!'
        ))
      }
    }

  } else {
    annotation_col <- NA
  }

  data[data > maximum] <- maximum
  if (minimum  != -Inf) {
    data[data < minimum] <-  minimum
  }

  cor_data = F
  dist_method = c('euclidean', "manhattan", "maximum", "canberra", "binary", "minkowski")

  if(scale == "row"){
    data_sd <- apply(data, 1, sd)
    data <- data[data_sd != 0, ]
  }

  if (correlation_plot  %in% c("row", "Row")) {
    if (clustering_distance_rows  == "pearson") {
      row_cor = cor(t(data))
    } else if (clustering_distance_rows  == "spearman") {
      row_cor = cor(t(data), method = "spearman")
    } else {
      if (clustering_distance_rows %in% dist_method){
        row_cor = as.data.frame(as.matrix(dist(data, method = clustering_distance_rows)))
      } else {
        row_cor = as.data.frame(as.matrix(vegan::vegdist(data, method = clustering_distance_rows)))
      }
    }
    data = round(row_cor, 3)
    annotation_col = annotation_row
    cor_data = T
  } else if (correlation_plot %in% c("col", "Column")) {
    # Do not know why add this!
    # Comment out
    # data_mad <- apply(data, 1, mad)
    # data <- data[data_mad > 0.1, ]
    if (clustering_distance_cols == "pearson") {
      col_cor = cor(data)
    } else if (clustering_distance_cols == "spearman") {
      col_cor = cor(data, method = "spearman")
    }  else {
      if (clustering_distance_cols %in% dist_method){
        col_cor = as.data.frame(as.matrix(dist(data, method = clustering_distance_rows)))
      } else {
        col_cor = as.data.frame(as.matrix(vegan::vegdist(data, method = clustering_distance_rows)))
      }
    }
    data = round(col_cor, 3)
    cor_data = T
    annotation_row = annotation_col
  }

  #print(data)
  # filter abnormal lines
  data_sd <- apply(data, 1, sd)
  if(any(data_sd==0)){
    stop("Wrong correlation method for this type of data. Please choose another method.")
  }



  if (width == 0 && height == 0) {
    height = nrow(data)
    width = ncol(data) * 1.1

    if (xtics_angle == 0) {
      width = width * 1.5
    }

  if (class(annotation_row) == "data.frame") {
    width = width + ncol(annotation_row)
    width = width * 1.1
  }

  if (class(annotation_col) == "data.frame") {
    height = height + ncol(annotation_col)
    width = width * 1.1
  }

  if (cluster_rows) {
    width = width + 4
  }

  if (cluster_cols) {
    height = height + 4
  }


    if (width < 8) {
      width = 8
    } else if (width < 20) {
      width = 8 + (width - 8) / 4
    } else if (width < 100) {
      width = 10 + (width - 20) / 5
    } else {
      width = 30
    }

    if (height < 10) {
      height = 8
    } else if (height < 20) {
      height = 8 + (height - 8) / 4
    } else if (height < 100) {
      height = 11 + (height - 20) / 5
    } else {
      height = 30
    }
  }






  if (sp.is.null(manual_annotation_colors_sidebar)) {
    manual_annotation_colors_sidebar = NA
  } else if (class(manual_annotation_colors_sidebar) == "character") {
    # Transfer string to R code
    manual_annotation_colors_sidebar = eval(parse(text = paste(
      "list(", manual_annotation_colors_sidebar, ")"
    )))
  }

  #print(manual_annotation_colors_sidebar)




  if (nrow(data) < 3) {
    cluster_rows = FALSE
    cluster_cols = FALSE
  }

  if (ncol(data) < 3) {
    cluster_cols = FALSE
    cluster_rows = FALSE
  }

  if (nrow(data) > 5000 & correlation_plot == "None") {
    cluster_rows = FALSE
  }

  if (ncol(data) > 5000 & correlation_plot == "None") {
    cluster_cols = FALSE
  }

  cluster_rows_results = cluster_rows
  cluster_cols_results = cluster_cols

  #if (height != 0) {
  #  height = height
  #}

  #if (width != 0) {
  #  width = width
  #}

  row_order = rownames(data)
  col_order = colnames(data)

  if (cluster_rows) {
    if (clustering_distance_rows == "pearson") {
      if (!cor_data) {
        row_cor = cor(t(data))
      } else {
        row_cor = data
      }
      row_dist <- as.dist(1 - row_cor)
      # Do not remember when this will happen
      if (any(is.na(row_cor))) {
        row_dist = dist(data)
      }
    } else if (clustering_distance_rows == "spearman") {
      if (!cor_data) {
        row_cor = cor(t(data), method = "spearman")
      } else {
        row_cor = data
      }
      row_dist <- as.dist(1 - row_cor)
      if (any(is.na(row_cor))) {
        row_dist = dist(data)
      }
    } else {
      if (!cor_data) {
        if (clustering_distance_rows %in% dist_method){
        row_dist = dist(data, method = clustering_distance_rows)
        } else {
          row_dist = vegan::vegdist(data, method = clustering_distance_rows)
        }
      } else {
        row_cor = data
        row_dist <- as.dist(1 - row_cor)
        if (any(is.na(row_cor))) {
          row_dist = dist(data)
        }
      }
    }
    cluster_rows_results = hclust(row_dist, method = clustering_method)

    if(sp.is.null(cluster_rows_variable)){
      sv = svd(data)$v[,1]
    } else {
      sv = annotation_row[[cluster_rows_variable]]
      if(remove_cluster_rows_variable_in_annorow){
        annotation_row[[cluster_rows_variable]] <- NULL
      }
      if(length(annotation_row)==0){
        annotation_row = NULL
      }
    }

    #print(sv)
    dend = reorder(as.dendrogram(cluster_rows_results), wts=sv)
    cluster_rows_results <- as.hclust(dend)
    row_order = cluster_rows_results$order
  }

  if (cluster_cols) {
    if (clustering_distance_cols == "pearson") {
      if (!cor_data) {
        col_cor = cor(data)
      } else {
        col_cor = data
      }
      col_dist <- as.dist(1 - col_cor)
      if (any(is.na(col_cor))) {
        col_dist = dist(t(data))
      }
    } else if (clustering_distance_cols  == "spearman") {
      if (!cor_data) {
        col_cor = cor(data, method = "spearman")
      } else {
        col_cor = data
      }
      col_dist <- as.dist(1 - col_cor)
      if (any(is.na(col_cor))) {
        col_dist = dist(t(data))
      }
    } else {
      if (!cor_data) {
        if (clustering_distance_cols %in% dist_method){
          col_dist = dist(t(data), method = clustering_distance_cols)
        } else {
          col_dist = vegan::vegdist(t(data), method = clustering_distance_cols)
        }
      } else {
        col_cor = data
        col_dist <- as.dist(1 - col_cor)
        if (any(is.na(col_cor))) {
          col_dist = dist(t(data))
        }
      }
    }
    cluster_cols_results = hclust(col_dist, method = clustering_method)
    if(sp.is.null(cluster_cols_variable)){
      sv = svd(data)$v[,1]
    } else {
      sv = annotation_col[[cluster_cols_variable]]

      if(remove_cluster_cols_variable_in_annocol){
        annotation_col[[cluster_cols_variable]] <- NULL
      }
      if(length(annotation_col) == 0){
        annotation_col = NULL
      }
    }

    dend = reorder(as.dendrogram(cluster_cols_results), wts=sv)
    cluster_cols_results <- as.hclust(dend)

    col_order = cluster_cols_results$order
  }


  if (correlation_plot!="None") {
    if (cluster_rows) {
      cluster_cols_results = cluster_rows_results
      col_order = row_order
    } else if (cluster_cols) {
      cluster_rows_results = cluster_cols_results
      row_order = col_order
    }
  }

  if (!is.na(cutree_rows)){
    if(mode(cluster_rows_results) != "logical"){
      data_row_cluster = as.data.frame(cutree(cluster_rows_results, cutree_rows))
      colnames(data_row_cluster) <- "Row_cluster"
      data_row_cluster$Row_cluster <- paste0("C", data_row_cluster$Row_cluster)
    }
  }

  if (anno_cutree_rows){
    if(is.na(annotation_row)){
      annotation_row = data_row_cluster
    } else {
      print(annotation_row)
      print(data_row_cluster)
      annotation_row = cbind(annotation_row, data_row_cluster)
    }
  }

  if (!is.na(cutree_cols)){
    if(mode(cluster_cols_results) != "logical"){
      data_col_cluster = as.data.frame(cutree(cluster_cols_results, cutree_cols))
      colnames(data_col_cluster) <- "Col_cluster"
      data_col_cluster$Col_cluster <- paste0("C", data_col_cluster$Col_cluster)
    }
  }

  if (anno_cutree_cols){
    if(is.na(annotation_col)){
      annotation_col = data_col_cluster
    } else {
      annotation_col = cbind(annotation_col, data_col_cluster)
    }
  }

  if(!is.na(filename)){

    data_order = data[row_order, col_order]
    sp_writeTable(data_order, file = paste0(filename, ".reordered.txt"))

    if (!is.na(cutree_rows)){
      if(mode(cluster_rows_results) != "logical"){
        data_row_cluster <- data_row_cluster[row_order, ,drop=F]
        sp_writeTable(
          data_row_cluster,
          file = paste0(filename, ".row_cluster.txt")
        )
      }
    }

    if (!is.na(cutree_cols)){
      if(mode(cluster_cols_results) != "logical"){
        data_col_cluster <- data_col_cluster[col_order, , drop=F]
        sp_writeTable(
          data_col_cluster,
          file = paste0(filename, ".col_cluster.txt")
        )
      }
    }

  }


  gt <- pheatmap::pheatmap(
    data,
    kmean_k = NA,
    color = manual_color_vector,
    scale = scale ,
    border_color = NA,
    cluster_rows = cluster_rows_results,
    cluster_cols = cluster_cols_results,
    cutree_rows = cutree_rows,
    cutree_cols = cutree_cols,
    kmeans_k = kclu,
    breaks = breaks,
    legend_breaks = legend_breaks,
    legend_labels = legend_labels,
    xtics_angle = xtics_angle,
    clustering_method = clustering_method ,
    clustering_distance_rows = clustering_distance_rows ,
    clustering_distance_cols = clustering_distance_cols ,
    show_rownames = ytics ,
    show_colnames = xtics ,
    main = title ,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = manual_annotation_colors_sidebar,
    fontsize = fontsize ,
    filename = filename,
    width = width,
    height = height,
    ...
  )

  if (saveppt){
  eoffice::topptx(gt, filename = paste0(filename,".pptx"))
  }
  gt
}
