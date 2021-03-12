#' A visualization of a correlation matrix.
#'
#' @param data Matrix or data file (with header line, the first column will be treated as row names, tab separated).
#' @param method The visualization method of correlation matrix to be used. Currently, it supports seven methods, named "circle" (default), "square", "ellipse", "number", "pie", "shade" and "color".
#' @param type Type, "full" (default), "upper" or "lower", display full matrix, lower triangular or upper triangular matrix.
#' @param bg The background color.
#' @param title Title of the graph.
#' @param is.corr Logical, whether the input matrix is a correlation matrix or not. We can visualize the non-correlation matrix by setting is.corr = FALSE.
#' @param diag Logical, whether display the correlation coefficients on the principal diagonal.
#' @param addCoef.col Color of coefficients added on the graph. If NULL (default), add no coefficients.
#' @param addCoefasPercent Logic, whether translate coefficients into percentage style for spacesaving.
#' @param order The ordering method of the correlation matrix. "original" for original order (default). "AOE" for the angular order of the eigenvectors. "FPC" for the first principal component order. "hclust" for the hierarchical clustering order. "alphabet" for alphabetical order.
#' @param hclust.method The agglomeration method to be used when order is hclust. This should be one of "ward", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param tl.cex Numeric, for the size of text label (variable names).
#' @param tl.col The color of text label.
#' @param tl.srt Numeric, for text label string rotation in degrees.
#' @param cl.pos Character or logical, position of color labels; If character, it must be one of "r" (default if type=="upper" or "full"), "b" (default if type=="lower") or "n", "n" means don't draw colorlabel.
#' @param cl.align.text "l", "c" (default) or "r", for number-label in colorlabel, "l" means left, "c" means center, and "r" means right.
#' @param cl.lim The limits (x1, x2) in the colorlabel.
#' @inheritParams base_plot_save
#' @param ... Other parameters given to base_plot_save
#' @param ...
#'
#' @return a grid object
#' @export
#'
#' @examples data<-matrix(rnorm(100),nrow=20)
#' rownames(data)<- paste0("corrtest", 1:20)
#' colnames(data)<- paste0("zhou", 1:5)
#' sp_corrplot(data, cl.align.text = "l")
#' sp_corrplot(data,method="pie",type="lower",title="corrplot", cl.align.text = "l", mar=c(0,0,1,0))
#' sp_corrplot(data,method="ellipse",type="upper", cl.align.text = "l", tl.pos="td")


sp_corrplot <-
  function(data,
           method = "circle",
           type = "full",
           # col_font = "black",
           bg = "white",
           title = NULL,
           is.corr = FALSE,
           diag = TRUE,
           addCoef.col=NULL,
           addCoefasPercent = FALSE,
           order = "original",
           hclust.method = NULL,
           tl.pos	= "n",
           tl.cex	= 1,
           tl.col	= "red",
           tl.srt=90,
           cl.pos	= NULL,
           cl.lim	=NULL,
           cl.align.text = "c",
           saveplot = NULL,
           mar = c(0, 0, 0, 0),
           ...) {

    if (class(data)[1] == "character") {
      data <- as.matrix(sp_readTable(data, row.names = 1))
    }

    if (is.corr == FALSE){
      rcorr_data <- rcorr(as.matrix(data))
      data <- rcorr_data$r
      correlation_test <- rcorr_data$P

      suppressWarnings(write.table(
        correlation_test,
        file = "correlation test.txt",
        sep = "\t",
        quote = F,
        row.names = TRUE
      ))
    }

    if (!sp.is.null(saveplot)) {
      base_plot_save(saveplot, ...)
    }

    corrplot::corrplot(
      data,
      method = method,
      type = type,
      # col = col_font,
      bg = bg,
      title = title,
      is.corr = is.corr,
      diag = diag,
      addCoef.col=addCoef.col,
      addCoefasPercent = addCoefasPercent,
      order = order,
      hclust.method = hclust.method,
      tl.pos	= tl.pos,
      tl.cex	= tl.cex,
      tl.col	= tl.col,
      tl.srt=tl.srt,
      cl.pos	= cl.pos,
      cl.lim	=cl.lim,
      cl.align.text=cl.align.text,
      mar = mar
    )
    if (!sp.is.null(saveplot)) {
      dev.off()
    }
  }
