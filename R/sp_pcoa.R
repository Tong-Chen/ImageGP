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
#' For `scale`: row scale.
#' For `sqrt` and `log2`, just as the words.
#' @param group_variable The variable for grouping points to generate normal data confidence ellipses. Optional.
#' @param draw_ellipse Default 'auto' means to enclose points in a polygon if one group with
#' less than 4 points. If there are more than 4 points for all groups, confidence ellipses would be draw.
#' Accept `confidence ellipse` to draw confidence ellipses for all conditions even though they would
#' not be draw.
#' Accept `no` to remove ellipse or other polygens.
#' @param label_variable The variable for text used to label points. Optional.
#' Specially supplying `Row.names` would label sample with their names.
#' @param check_significance Check if the centroids and dispersion of the groups
#' as defined by measure space are equivalent for all groups.
#' @param coord_fixed  When True (the default) ensures that one unit on the x-axis is
#' the same length as one unit on the y-axis.
#' @param check_paired_significance Paired-check for each two groups.
#' @inheritParams sp_boxplot
#' @inheritParams stats::cmdscale
#' @inheritParams ggplot2::stat_ellipse
#' @inheritParams dataFilter2
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
                    top_n = 1,
                    statistical_value_type = mad,
                    binary_dissimilarity_index = F,
                    data_transform = 'auto',
                    group_variable = NULL,
                    color_variable = NULL,
                    color_variable_order = NULL,
                    shape_variable = NULL,
                    shape_variable_order = NULL,
                    size_variable = NULL,
                    size_variable_order = NULL,
                    label_variable = NULL,
                    label_variable_order = NULL,
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
                    coord_fixed = T,
                    ...) {
  if ("character" %in% class(data)) {
    file <- data
    data <- sp_readTable(data, row.names = 1)
    if (input_type == "normalized_OTUtable" &&
        (!'dist' %in% class(data))){
      data <- dataFilter2(data, top_n = top_n,
                          statistical_value_type = statistical_value_type)
    }

  } else if ('data.frame' %in% class(data) |
             'dist' %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  if (class(metadata) == "character") {
    metadata <- sp_readTable(metadata, row.names = 1)
  } else if ('data.frame' %in% class(metadata) |
             'dist' %in% class(metadata)) {
    stop("Unknown input format for `metadata` parameter.")
  }

  # Keep same columns of data with rows of metadata
  matchedL <- match_two_df(data, metadata, way = "col-row")

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

  distance_algorithm = ""
  if (input_type == "normalized_OTUtable" &&
      (!'dist' %in% class(data))) {
    data <- t(data)
    # copy and modified from vegan::metaMDSdist
    if (data_transform == "auto") {
      xam <- max(data, na.rm=F)
      if (xam > 50) {
        data <- sqrt(data)
      }
      if (xam > 9) {
        data <- vegan::wisconsin(data)
      }
    } else if (data_transform == "sqrt") {
        data <- sqrt(data)
    } else if (data_transform == "log2") {
      log_add = sp_determine_log_add(data)
      data <- data + log_add
      data <- log2(data)
    }  else if (data_transform == "scale") {
      data <- data[var(data) != 0, ]
      data_scale <- t(apply(data, 1, scale))
      rownams(data_scale) <- rownames(data)
      colnams(data_scale) <- colnames(data)
      data <- data_scale
    } else if (data_transform != "None") {
      data <- vegan::decostand(data, method = data_transform)
    }

    dist_matrix <- vegan::vegdist(data, method = dissimilarity_index,
                                  binary = binary_dissimilarity_index,
                                  na.rm=TRUE)
    distance_algorithm <- dissimilarity_index
  } else {
    dist_matrix <- as.dist(data)
  }


  # 设置维数
  ndimensions = k
  num_samples <- nrow(data)
  if (ndimensions >= num_samples) {
    ndimensions <- num_samples - 1
  }

  # 计算pcoa
  pcoa <- cmdscale(dist_matrix, k = ndimensions, eig = T)



  # 整理pcoa结果
  pcoa_points <- as.data.frame(pcoa$points)
  sum_eig <- sum(pcoa$eig)
  eig_percent <- round(pcoa$eig / sum_eig * 100, 1)

  colnames(pcoa_points) <- paste0("PCoA", 1:ndimensions)
  print(pcoa_points)

  data <- cbind(pcoa_points, metadata)

  data$Row.names <- row.names(metadata)
  rownames(data) <- data$Row.names

  data_colnames <- colnames(data)
  #print(data_colnames)
  #print(data)

  if (!sp.is.null(color_variable)) {
    if (sp.is.null(group_variable)) {
      group_variable =  color_variable
    }
    if (!(color_variable %in% data_colnames)) {
      stop(paste(color_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, color_variable, color_variable_order)
  }


  if (!sp.is.null(shape_variable)) {
    if (sp.is.null(group_variable)) {
      group_variable =  shape_variable
    }
    if (!(shape_variable %in% data_colnames)) {
      stop(paste(shape_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, shape_variable, shape_variable_order)
    data[[shape_variable]] <- as.factor(data[[shape_variable]])
    shapes <- generate_shapes(data, shape_variable)
  }

  if (!sp.is.null(size_variable)) {
    if (!(size_variable %in% data_colnames)) {
      stop(paste(size_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, size_variable, size_variable_order)
  }



  if (!sp.is.null(group_variable) &&
      group_variable != color_variable &&
      group_variable != shape_variable) {
    if (!(group_variable %in% data_colnames)) {
      stop(paste(group_variable, 'must be column names of data!'))
    }
    #data = sp_set_factor_order(data, group_variable, group_variable_order)
  }

  if (draw_ellipse == 'auto') {
    if (all(table(data[[group_variable]]) > 4)) {
      draw_ellipse = "confiden ellipse"
    } else {
      # library(ggalt)
      draw_ellipse = "encircle"
    }
  }

  library(ggplot2)

  group_variable_en = sym(group_variable)

  if (!sp.is.null(filename)) {
    sp_writeTable(
      data,
      file = paste0(filename, ".pcoas.txt"),
      keep_rownames = F
    )
  }

  analysis_label = "PCoA"
  if (input_type == "normalized_OTUtable"){
    if (distance_algorithm != ""){
      analysis_label = paste("PCoA", "(", distance_algorithm, ")", sep="")
    }

    if (dissimilarity_index == "euclidean"){
      analysis_label = "PCA"
    }
  }


  p <- ggplot(data, aes(x = PCoA1, y = PCoA2, group = !!group_variable_en)) +
    labs(
      x = paste(analysis_label, "1 (", eig_percent[1], "%)", sep = ""),
      y = paste(analysis_label, "2 (", eig_percent[2], "%)", sep = ""),
      title = title
    ) + geom_point()

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

  if (!sp.is.null(label_variable)) {
	# For cloud platform usages.
    if (!(label_variable %in% data_colnames)) {
      label_variable = "Row.names"
    }
    label_variable_en = sym(label_variable)
    library(ggrepel)
    p <-
      p + geom_text_repel(
        aes(label = !!label_variable_en),
        show.legend = F,
        max.overlaps = 100
      )
  }

  if (draw_ellipse == "encircle") {
    p <-
      p + geom_encircle(alpha = 0.2,
                        show.legend = F,
                        aes(fill = !!group_variable_en))
    p <-
      sp_manual_fill_ggplot2(p, data, group_variable, manual_color_vector)

  } else if (draw_ellipse == "confiden ellipse") {
    p <-
      p + stat_ellipse(
        level = level,
        type = type,
        na.rm = TRUE
      )
  }

  if (check_significance) {
    pcoa_adonis2 <-
      adonis2(as.formula(paste("dist_matrix", "~", group_variable)),
              data = metadata,
              permutations = 5999)

    dispersion <-
      betadisper(dist_matrix, group = metadata[[group_variable]])
    dispersion_test <- permutest(dispersion)
    dispersion_test_p <- dispersion_test$tab$`Pr(>F)`[1]

    title <- paste0(
      "adonis R2: ",
      round(pcoa_adonis2$R2, 2),
      "\nadonis P-value: ",
      round(pcoa_adonis2$`Pr(>F)`,6),
      "\ndispersion P-value: ",
      round(dispersion_test_p,6)
    )

    if (check_paired_significance) {
      # devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
      library(pairwiseAdonis)
      pairwise.adonis <-
        pairwise.adonis(
          x = dist_matrix,
          factors = metadata[[group_variable]],
          p.adjust.m = "BH",
          reduce = NULL,
          perm = 5999
        )
      if (!sp.is.null(filename)) {
        sp_writeTable(
          pairwise.adonis,
          file = paste0(filename, ".pairwiseAdonis.txt"),
          keep_rownames = F
        )
      }

      tukeyHSD <- TukeyHSD(dispersion)$group

      if (!sp.is.null(filename)) {
        sp_writeTable(
          tukeyHSD,
          file = paste0(filename, ".pairwiseDispersionCheck.txt")
        )
      }

    }
  }



  p <- sp_ggplot_layout(
    p,
    legend.position = legend.position,
    extra_ggplot2_cmd = extra_ggplot2_cmd,
    filename = filename,
    title = title,
    ...
  )

  if (coord_fixed){
    p <- p +  coord_fixed(1)
  }

  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }
  p
}
