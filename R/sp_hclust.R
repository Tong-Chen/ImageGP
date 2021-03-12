#' Hierarchical cluster diagram
#'
#' @param data A matrix file or an object.
#' @param method Clustering method :"ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"
#' @param thresholdZ.k Threshold for defining outliers. First compute the overall
#' corelation of one sample to other samples. Then do Z-score transfer for all
#' correlation values. The samples with corelation values less than given value
#' would be treated as outliers.
#' Default -2.5 meaning -2.5 std.
#' @param ...
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' x = runif(10)
#' y = runif(10)
#' data=cbind(x, y)
#' rownames(data) = paste("exam", 1:10)
#' sp_hclust(data)
#'
sp_hclust <- function (data,
                       method = "average",
                       thresholdZ.k = -2.5,
                       saveplot = NULL,
                       debug = FALSE,
                       ...) {
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }
  if (class(data) == "character") {
  datExpr <- sp_readTable(data, row.names = NULL)
  } else {
    datExpr <- data
  }
  A = WGCNA::adjacency(t(datExpr), type = "distance")
  # this calculates the whole network connectivity
  k = as.numeric(apply(A, 2, sum)) - 1
  # standardized connectivity
  Z.k = scale(k)
  # Designate samples as outlying if their Z.k value is below the threshold
  # thresholdZ.k = -5  # often -2.5

  if (thresholdZ.k > 0) {
    cat("\tThe program will transfer positive thresholdZ.k to their negative values.\n")
    thresholdZ.k = -1 * thresholdZ.k
  }

  cat("\tThreshold for detecting outlier samples are",
      thresholdZ.k,
      "\n")
  # the color vector indicates outlyingness (red)
  outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")

  # calculate the cluster tree using flahsClust or hclust
  sampleTree = hclust(as.dist(1 - A), method = method)

  if (!sp.is.null(saveplot)) {
    base_plot_save(saveplot, ...)
  }
  plotDendroAndColors(
    sampleTree,
    groupLabels = names(outlierColor),
    colors = outlierColor,
    main = "Sample dendrogram"
  )
  if (!sp.is.null(saveplot)) {
    dev.off()
  }

}
