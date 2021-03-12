#' Hierarchical cluster diagram
#'
#' @param data A data frame.
#' @param group_variable Specifies a column as group.
#' @param branch_order Specify branch order.
#' @param method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param k Number of groups.
#' @param labels_size Labels size.
#' @param shape Tree or circles.
#' @param pic_title Title of the graph.
#' @param pic_flip TRUE for horizontal, FALSE for vertical. Default TRUE.
#' @param node_size Node size.
#' @param legend_site Legend site. Optional, "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
#' @param ...
#'
#' @return
#' @export
#'
#' @examples data<-matrix(rnorm(200),nrow = 50)
#' rownames(data)<- paste0("dendextend", 1:50)
#' colnames(data)<- paste0("zhou", 1:4)
#' sp_dendextend(data = data,k = 3,labels_size = 0.3)
#'
#'
#' data<- data.frame(ID = letters[1:5], apple = runif(50), banana = runif(50), watermelon = runif(50))
#' sp_dendextend(data = data, k = 5, method = "single", shape = "circle")
#'


sp_dendextend <- function(data,
                          group_variable = NULL,
                          branch_order = NULL,
                          method = "complete",
                          k = k,
                          labels_size = 0.5,
                          shape = "tree",
                          pic_title = NULL,
                          pic_flip = TRUE,
                          node_size = 0.007,
                          legend_site = "topleft",
                          saveplot = NULL,
                          ...) {
  if (class(data)[1] == "character") {
    data <- sp_readTable(data, row.names = NULL)
  }

  if (!sp.is.null(group_variable)) {
    # subset(data, select =-"Species")
    matrix_data <-
      data[, -which(names(data) %in% c(group_variable))]
    group_data <- data[, group_variable]
  }

  dist_data <- dist(matrix_data)
  hc_data <- hclust(dist_data, method = method)

  dend <- as.dendrogram(hc_data)

  # order it the closest we can to the order of the observations:
  if (!sp.is.null(branch_order)) {
    dend <-
      dendextend::rotate(dend, order = as.character(branch_order))
  }

  # Color the branches based on the clusters:
  dend <-
    dendextend::color_branches(dend, k = k)#, groupLabels=iris_species)

  if (!sp.is.null(group_variable)) {
    # Manually match the labels, as much as possible, to the real classification of the flowers:
    labels_colors(dend) <-
      rainbow_hcl(k)[sort_levels_values(as.numeric(as.factor(data[, group_variable]))[order.dendrogram(dend)])]
  }

  if (!sp.is.null(group_variable)) {
    # We shall add the flower type to the labels:
    labels(dend) <-
      paste(as.character(group_data)[order.dendrogram(dend)],
            "(", labels(dend), ")",
            sep = "")
  }

  dend <- hang.dendrogram(dend, hang_height = 0.1)
  dend <- dendextend::set(dend, "labels_cex", labels_size)

  if (!sp.is.null(saveplot)) {
    base_plot_save(saveplot, ...)
  }
  if (shape == "tree") {
    plot(
      dend,
      main = pic_title,
      horiz = pic_flip,
      nodePar = list(cex = node_size)
    )
    if (!sp.is.null(group_variable) || !sp.is.null(legend_site)) {
      group_data <- rev(levels(as.factor(data[, group_variable])))
      legend(legend_site, legend = group_data, fill = rainbow_hcl(k))
    }
  } else {
    par(mar = rep(0, 4))
    circlize_dendrogram(dend)
  }
  if (!sp.is.null(saveplot)) {
    dev.off()
  }
}
