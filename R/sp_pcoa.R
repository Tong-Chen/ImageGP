#' Generating pcoa plot
#'
#' @param data Data file (or data.frame) With the first row as the header line and the first column as the row-name.
#' Columns are separated by one `tab`. Each row represents one variable (normally genes, OTU).
#' Each column represents one sample. The numbers represent gene expression abundance or OTU abundance or
#' other abundances, and should be normalized.
#' @param metadata Metadata file (or data.frame) with sample attributes like group information.
#' The first column is the same as the first row of value given to parameter `data`.
#' These attributes would be used as `color`, `size`, `shape` variables in the plot.
#' If not supplied, each sample will be treated as one group.
#' @param dissimilarity_index Dissimilarity index, partial match to "manhattan", "euclidean",
#' "canberra", "clark", "bray" (default, meaning Bray–Curtis), "kulczynski", "jaccard", "gower", "altGower",
#' "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao",
#' "mahalanobis", "chisq" or "chord".
#' Gower, Bray–Curtis, Jaccard and Kulczynski indices are good in detecting underlying
#' ecological gradients (Faith et al. 1987).
#' Morisita, Horn–Morisita, Binomial, Cao and Chao indices should be able to handle different
#' sample sizes (Wolda 1981, Krebs 1999, Anderson & Millar 2004),
#' and Mountford (1962) and Raup-Crick indices for presence–absence data should be able
#' to handle unknown (and variable) sample sizes.
#' @param binary_dissimilarity_index Perform presence/absence standardization
#' before computing `dissimilarity_index`. Default `FALSE`, accept `TRUE`.
#' @param input_type The input data is OTU table (`normalized_OTUtable`) or
#' a distance matrix (`distance_matrix`).
#' @param data_transform Methods for transforming data. Default 'auto'. Accept 'None',
#' For `auto`: If the data values are larger than common abundance class scales
#' (here is `9`), the function performs a Wisconsin double standardization (\code{\link[vegan]{wisconsin}}).
#' If the values look very large, the function also performs \code{\link{sqrt}} transformation.
#' For `None`: No transformation would be performed.
#' For `total`: Compute relative abundance of OTUs/Genes in each sample.
#' For `hellinger`: square root of `method = "total"`.
#' @param group_variable The variable for grouping points to generate normal data confidence ellipses. Optional.
#' @param draw_ellipse Default 'auto' means to enclose points in a polygon if one group with
#' less than 4 points. If there are more than 4 points for all groups, confidence ellipses would be draw.
#' Accept `confidence ellipse` to draw confidence ellipses for all conditions even though they would
#' not be draw.
#' @param label_variable The variable for text used to label points. Optional.
#' @param check_significance Check if the centroids and dispersion of the groups
#' as defined by measure space are equivalent for all groups.
#' @param check_paired_significance Paired-check for each two groups.
#' @inheritParams sp_boxplot
#' @inheritParams stats::cmdscale
#' @inheritParams ggplot2::stat_ellipse
#' @inheritParams sp_scatterplot
#' @param ... Parameters given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' ## End(Not run)
#'
sp_pcoa <- function(data,
                    metadata,
                    input_type = "normalized_OTUtable",
                    dissimilarity_index = "bray",
                    k = 3,
                    binary_dissimilarity_index = F,
                    data_transform = 'auto',
                    group_variable = NULL,
                    color_variable = NULL,
                    color_variable_order = NULL,
                    shape_variable = NULL,
                    shape_variable_order = NULL,
                    size_variable = NULL,
                    size_variable_order = NULL,
                    label_variable = FALSE,
                    legend.position = 'right',
                    draw_ellipse = 'auto',
                    manual_color_vector = NULL,
                    title = NULL,
                    label_font_size = NULL,
                    debug = FALSE,
                    type = 't',
                    level = 0.95,
                    filename = NULL,
                    extra_ggplot2_cmd = NULL,
                    check_significance = T,
                    check_paired_significance = T,
                    ...) {


  if (class(data) == "character") {
    file <- data
    data <- sp_readTable(data, row.names = 1)
  } else if ('data.frame' %in% class(data) | 'dist' %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  if (class(metadata) == "character") {
    metadata <- sp_readTable(metadata, row.names = 1)
  } else if ('data.frame' %in% class(metadata) | 'dist' %in% class(metadata)) {
    stop("Unknown input format for `metadata` parameter.")
  }

  # Keep same columns of data with rows of metadata
  matchedL <- match_two_df(data, metadata, way="col-row")

  data = matchedL$df1
  metadata = matchedL$df2

  # metadata_colnames <- colnames(metadata)
  # if(!sp.is.null(group_variable) && (!group_variable %in% metadata_colnames) ){
  #   stop(paste(group_variable, 'must be column names of data!'))
  # }
  #
  # if(!sp.is.null(color_variable) && (!color_variable %in% metadata_colnames) ){
  #   stop(paste(color_variable, 'must be column names of data!'))
  # }

  if(input_type == "normalized_OTUtable" && (!'dist' %in% class(data)) ){
    data <- t(data)
    # copy and modified from vegan::metaMDSdist
    if(data_transform == "auto"){
      xam <- max(data)
      if (xam > 50) {
        data <- sqrt(data)
      }
      if (xam > 9) {
        data <- vegan::wisconsin(data)
      }
    } else if(data_transform != "None"){
      data <- vegan::decostand(data, method=data_transform)
    }

    dist_matrix <- vegan::vegdist(data, method=dissimilarity_index,
                           binary=binary_dissimilarity_index)
  } else {
    dist_matrix <- as.dist(data)
  }


  # 设置维数
  ndimensions = k
  num_samples <- nrow(data)
  if(ndimensions >= num_samples){
    ndimensions <- num_samples - 1
  }

  # 计算pcoa
  pcoa <- cmdscale(dist_matrix, k=ndimesions, eig=T)



  # 整理pcoa结果
  pcoa_points <- as.data.frame(pcoa$points)
  sum_eig <- sum(pcoa$eig)
  eig_percent <- round(pcoa$eig/sum_eig*100,1)

  colnames(pcoa_points) <- paste0("PCoA", 1:ndimensions)

  data <- cbind(dune_pcoa_points, metdata)

  data_colnames <- colnames(data)



  # if (!is.numeric(data[[yvariable]]) &&
  #     !sp.is.null(yvariable_order)) {
  #   data = sp_set_factor_order(data, yvariable, yvariable_order)
  # }
  #
  # if (!is.numeric(data[[yvariable]]) &&
  #     !sp.is.null(xvariable_order)) {
  #   data = sp_set_factor_order(data, xvariable, xvariable_order)
  # }


  if (!sp.is.null(color_variable)) {
    if(sp.is.null(group_variable)){
      group_variable =  color_variable
    }
    if (!(color_variable %in% data_colnames)) {
      stop(paste(color_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, color_variable, color_variable_order)
  }


  if (!sp.is.null(shape_variable)) {
    if(sp.is.null(group_variable)){
      group_variable =  shape_variable
    }
    if (!(shape_variable %in% data_colnames)) {
      stop(paste(shape_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, shape_variable, shape_variable_order)
    shapes <- generate_shapes(data, shape_variable)
  }

  if (!sp.is.null(size_varaible)) {
    if (!(size_varaible %in% data_colnames)) {
      stop(paste(size_varaible, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, size_varaible, size_varaible_order)
  }

  if (!sp.is.null(group_variable) && group_variable != color_variable &&
      group_variable != shape_variable) {
    if (!(group_variable %in% data_colnames)) {
      stop(paste(group_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, group_variable, group_variable_order)
  }

  if(draw_ellipse == 'auto'){
    if(all(table(data[[group_variable]])>4)){
      draw_ellipse = "confiden ellipse"
    } else {
      # library(ggalt)
      draw_ellipse = "encircle"
    }
  }

  library(ggplot2)

  p <- ggplot(data, aes(x=PCoA1, y=PCoA2, group=group_variable)) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title=title) + geom_point()

  if (!sp.is.null(color_variable)) {
      color_variable_en = sym(color_variable)
      p <- p + aes(color = !!color_variable_en)
      p <-
        sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)
    }

  if (!sp.is.null(shape_variable)) {
    shape_variable_en = sym(shape_variable)
    p <- p + aes(shape = !!shape_variable_en) +
      scale_shape_manual(values = shapes)
  }

  if (!sp.is.null(size_variable)) {
    size_variable_en = sym(size_variable)
    p <- p + aes(size = !!size_variable_en)
  }

  if (!sp.is.null(label_variable)){
    label_variable_en = sym(label_variable)
    library(ggrepel)
    p <- p + geom_text_repel(aes(label=!!label_variable_en), show.legend = F,
                             max.overlaps = 100)
  }

  if(draw_ellipse == "encircle"){
    p <- p + geom_encircle(aes(fill=group_variable), alpha = 0.1, show.legend = F)
  } else{
    p <- p + stat_ellipse(level = level, type=type, na.rm = TRUE)
  }

  if(check_significance){
    pcoa_adonis2 <- adonis2(as.formula(paste("dist_matrix", "~", group_variable)),
                                             data = metadata,
                                             permutations = 5999)

    dispersion <- betadisper(dist_matrix, group=metadata[[group_variable]])
    dispersion_test <- permutest(dispersion)
    dispersion_test_p <- dispersion_test$tab$`Pr(>F)`[1]

    title <- paste0("adonis R2: ",round(pcoa_adonis2$R2,2),
                    "; adonis P-value: ", pcoa_adonis2$`Pr(>F)`,
                    "; dispersion P-value: ", dispersion_test_p)

    if(check_paired_significance){
      # devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
      library(pairwiseAdonis)
      dune.pairwise.adonis <- pairwise.adonis(x=dist_matrix, factors=metadata[[group_variable]],
                                              p.adjust.m = "BH",
                                              reduce = NULL,
                                              perm = 5999)
      if(!is.na(filename)){
        write.table(
          dune.pairwise.adonis,
          file = paste0(filename, ".pairwiseAdonis.txt"),
          sep = "\t",
          quote = F,
          col.names = T,
          row.names = F
        )
      }

    }
  }


  if (!sp.is.null(facet_variable)) {
    if (facet_singlecell_style) {
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
    coordinate_flip = coordinate_flip,
    filename = filename,
    additional_theme = additional_theme,
    ...
  )

  p <- p +  coord_fixed(1)




  p

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
