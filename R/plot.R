# www.ehbio.com/Training
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'


#' Boxplot for wide data frame (normal gene expression table or OTU abundance table)
#' used for estimating the overall distrbution of data.
#'
#' @param widedataframe A dataframe containing gene expression or OTU abundance with format
#' like generated by \code{\link{generateAbundanceDF}}.
#' @param saveplot Save plot to given file like "a.pdf", "b.png".
#' @param ylab Y axis title.
#' @param ... Other parameters given to \code{\link[ggplot2]{ggsave}}.
#'
#' @return A ggplot2 object.
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' widedataframe2boxplot(df)
#'
#' widedataframe2boxplot(df, saveplot="a.pdf", width=10, height=10, units=c("cm"))
#'
widedataframe2boxplot <- function(widedataframe, saveplot=NULL, ylab="", ...) {
  #widedataframe2 <- widedataframe
  widedataframe$id <- rownames(widedataframe)
  rlog_mat_melt <- reshape2::melt(widedataframe, id.vars = c('id'))
  p <- ggplot(rlog_mat_melt, aes(x=variable, y=value)) + geom_boxplot(aes(color=variable)) +
    geom_violin(aes(fill=variable), alpha=0.5) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(), legend.position = "none") + ylab(ylab)

  if(!is.null(saveplot)){
    ggsave(filename=saveplot, p, ...)
  }
  return(p)
}


#' Rankplot for given column.
#'
#' @param data A dataframe with effective row names.
#' @param order_col Specify which column would be used for plot. Default "log2FoldChange".
#' @param midpoint Specify the midpoint for color show. Default 0.
#' @param saveplot Save plot to given file "a.pdf", "b.png".
#' @param label Label points. Accept a number to get top x points to label (both directions).
#' Or a vector matched with rownames of dataframe to label specified points.
#' Or a named vector to label specified points with new names.
#' @param colorvector A vector with length 2 to speicfy low and high colors.
#' Or a vector with length 3 to specify low, middle and high colors.
#' Default \code{c("green","yellow","red")}.
#' @param alpha Transparency value.
#' @param width Picture width in "cm".
#' @param height Picture height in "cm".
#' @param ... Additional parameters given to \code{\link[pheatmap]{pheatmap}}
#' like "width", "height", .etc.
#'
#' @import ggplot2
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' a <- data.frame(log2FoldChange=rnorm(1000), row.names=paste0("ImageGP",1:1000))
#'
#' ## Raw plot
#'
#' rankPlot(a)
#'
#' ## Label top 10
#'
#' rankPlot(a, label=10)
#'
#' ## Label specified points
#'
#' rankPlot(a, label=c("ImageGP1","ImageGP2","ImageGP10"))
#'
#' ## Label specified points with new names
#'
#' b <- c("A","B","C")
#' names(b) <- c("ImageGP1","ImageGP2","ImageGP10")
#' rankPlot(a, label=b)
#'

rankPlot <- function(data, order_col="log2FoldChange",
                     midpoint=0,
                     saveplot=NULL,
                     # A number to get top x labels
                     # Or a vector of ids (should match rownames) to show
                     # Or a named vector to label new symbols
                     label=NULL,
                     alpha=0.5,
                     colorvector=c("green","yellow","red"),
                     width=13.5, height=15, ...){
  data_line <- data[order(data[[order_col]]),,drop=F]
  row_num <- nrow(data_line)
  data_line$x <- 1:row_num
  data_line$xend <- 2:(row_num+1)
  if(!is.null(label)){
    if (length(label) == 1 & is.numeric(label)) {
      # Label top
      data_line$id__ct <- rownames(data_line)
      data_line[(label+1):(row_num-label),"id__ct"] <- NA
    } else {
      if(is.null(names(label))){
        # Label given ids
        data_line$id__ct <- NA
        data_line$id__ct <- label[match(rownames(data_line), label)]
      } else {
        # Label substitute ids
        data_line$id__ct <- NA
        data_line$id__ct <- label[match(rownames(data_line), names(label))]
      }
    }
  }


  p <- ggplot(data_line) +
    geom_segment(aes_string(x="x", xend="x",y=order_col,yend=0, color=order_col), alpha=alpha)
  #geom_linerange(aes_string(x="x", ymin=order_col,ymax=0, color=order_col), alpha=alpha)
  #geom_rect(aes_string(xmin="x", xmax="xend",ymin=order_col,ymax=0, fill=order_col), alpha=alpha)

  if(length(colorvector)==2){
    p <- p + scale_color_gradient(low=colorvector[1], high=colorvector[2])
  } else if(length(colorvector)==3){
    #names(colorvector) <- c("low","mid","high")
    #colorvector <- as.list(colorvector)
    #p <- p + scale_color_gradient2(colorvector, midpoint = midpoint)
    p <- p + scale_color_gradient2(low=colorvector[1], mid=colorvector[2],
                                   high=colorvector[3], midpoint = midpoint)
  } else {
    stop("Illegel color vector.")
  }
  p <- p + theme_classic() + geom_hline(yintercept = midpoint, linetype="dotted")

  if(!is.null(label)){
    # 标记基因名字
    p <- p + ggrepel::geom_text_repel(aes_string(x="x", y=order_col, label="id__ct"))
  }

  if(!is.null(saveplot)){
    ggsave(plot=p, filename=saveplot, units=c("cm"),...)
  }

  return(p)

}


#' Volcano plot
#'
#' @param data A dataframe.
#' @param log2FoldChange Specify the columns containing log2 fold change.
#' Default "log2FoldChange" (suitable for DESeq2 result)
#' @param padj Specify the columns containing adjusted p-value.
#' Default "padj"  (suitable for DESeq2 result)
#' @param colour Specify colour variable. Normally the columns containing
#' labels to indicate if the genes are up-regulated or down-regulated or
#' no significant difference.
#' @param saveplot Save plot to given file "a.pdf", "b.png".
#' @param size A number of one column name to specify point size.
#' @param padjLimit Max allowed negative log10 transformed padj.
#' Default 10.
#' @param width Picture width in "cm".
#' @param height Picture height in "cm".
#' @param ... Additional parameters given to \code{\link[ggplot2]{ggsave}}.
#'
#' @import ggplot2
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#'
#' ### Generate test data
#' res_output <- data.frame(log2FoldChange=rnorm(3000), row.names=paste0("ImageGP",1:3000))
#' res_output$padj <- 20 ^ (-1*(res_output$log2FoldChange^2))
#' padj = 0.05
#' log2FC = 1
#' res_output$level <- ifelse(res_output$padj<=padj,
#'                            ifelse(res_output$log2FoldChange>=log2FC,
#'                                   paste("groupA","UP"),
#'                                   ifelse(res_output$log2FoldChange<=(-1)*(log2FC),
#'                                          paste("groupB","UP"), "NoDiff")) , "NoDiff")
#' head(res_output)
#'
#' volcanoPlot(res_output, colour="level")
#'
volcanoPlot <- function(data, log2FoldChange="log2FoldChange", padj="padj",
                        colour='red',
                        saveplot=NULL, size=1, padjLimit=10,
                        width=13.5, height=15, ...){
  data[[padj]] <- (-1) * log10(data[[padj]])

  #if (is.null(padjLimit)) {
  #  #padjLimit <- quantile(data[[padj]], probs=seq(0,1,0.1))[10]
  #  print(paste("Max allowed negative log10 padj", padjLimit))
  #}

  data[[padj]] <- replace(data[[padj]], data[[padj]]>padjLimit, padjLimit*1.001)

  boundary = ceiling(max(abs(data[[log2FoldChange]])))

  if (is.numeric(size)){
    show.legend = F
  }

  # p = ggplot(data, aes(x=!!ensym(log2FoldChange),y=!!ensym(padj),
  #  colour=!!ensym(level)))
  p = ggplot(data, aes_string(x=log2FoldChange,y=padj,colour=colour))
  if (is.numeric(size)){
    p <- p + geom_point(size=size, alpha=0.5)
  } else {
    p <- p + geom_point(aes_string(size=size), alpha=0.5)
  }

  p <- p +
    theme_classic() +
    xlab("Log2 transformed fold change") +
    ylab(paste("Negative Log10 transformed", padj)) +
    xlim(-1 * boundary, boundary) +
    theme(legend.position="top", legend.title=element_blank())

  if(!is.null(saveplot)){
    ggsave(plot=p, filename=saveplot, units=c("cm"),...)
  }
  return(p)

}




#' Compute and plot column correlation matrix. Normally used
#' to do sample corealtion of gene expression or OTU abundance matrix.
#'
#' @inheritParams Matrix2colCorrelation
#' @param saveplot Save plot to given file "a.pdf", "b.png".
#' @param width Picture width in "cm".
#' @param height Picture height in "cm".
#' @param ... Additional parameters given to \code{\link[ggplot2]{ggsave}}.
#'
#' @import ggplot2
#' @import grDevices
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' clusterSampleUpperTriPlot(df)
#'
clusterSampleUpperTriPlot <- function (mat, method="pearson", digits=4,
                                       cor_file=NULL, saveplot=NULL,
                                       width=13.5, height=15, ...){
  print("Performing sample clustering")

  pearson_cor_hc <- Matrix2colCorrelation(mat, method, digits, cor_file)
  pearson_cor <-  pearson_cor_hc$pearson_cor

  upper_tri <- get_upper_tri(pearson_cor)
  # Melt the correlation matrix
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

  col = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")))(100)

  # Create a ggheatmap
  p <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colours=col, name=paste(method,"correlation")) +
    theme_classic() +
    coord_fixed() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.justification = c(1, 0),
      legend.position = "top",
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 9, barheight = 1,
                                 title.position = "left"))

  if(!is.null(saveplot)){
    ggsave(plot=p, filename=saveplot, units=c("cm"),...)
  }
  return(p)
}

#' Generate suitable output graphics device by file suffix.
#'
#' @param saveplot Save plot to given file "a.pdf", "b.png".
#' @param ... Additional parameters given to plot output (\code{\link{pdf}}, \code{\link{png}},...) like "width", "height", .etc.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' base_plot_save("a.pdf")
#' # will simplify run (pdf("a.pdf))
#'
base_plot_save <- function(saveplot, ...) {
  # From https://github.com/raivokolde/pheatmap/blob/master/R/pheatmap.r
  # Get file type
  r = regexpr("\\.[a-zA-Z]*$", saveplot)
  if(r == -1) stop("Improper filename")
  ending = substr(saveplot, r + 1, r + attr(r, "match.length"))

  f = switch(ending,
             pdf = function(x, ...) pdf(x, paper="special", useDingbats=F, ...),
             png = function(x, ...) png(x, units = "in", res = 300, ...),
             jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
             jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
             tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
             bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
             stop("File type should be: pdf, png, bmp, jpg, tiff")
  )
  f(saveplot, ...)

}

#' Compute and plot column correlation matrix. Normally used
#' to do sample corealtion of gene expression or OTU abundance matrix.
#'
#' @inheritParams Matrix2colCorrelation
#' @inheritParams base_plot_save
#' @inheritParams dataFilter
#'
#'
#' @return Nothing
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' clusterSampleHeatmap2(df)
#'
clusterSampleHeatmap2 <- function (mat, method="pearson", digits=4,
                                   cor_file=NULL, saveplot=NULL, ...){
  print("Performing sample clustering")

  hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(100)

  mat <- dataFilter(mat, ...)

  pearson_cor_hc <- Matrix2colCorrelation(mat, method, digits, cor_file)
  pearson_cor <-  pearson_cor_hc$pearson_cor
  hc <-  pearson_cor_hc$hc

  if(!is.null(saveplot)) {
    base_plot_save(saveplot, ...)
  }

  gplots::heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
                    col=hmcol, margins=c(11,11), key=T,
                    main="Correlation plot")
  if(!is.null(saveplot)) {
    dev.off()
  }
}


#' Compute and plot column correlation matrix. Normally used
#' to do sample corealtion of gene expression or OTU abundance matrix.
#'
#' @inheritParams Matrix2colCorrelation
#' @param sampleAnno Add sample attribute to plot. Accept a dataframe with
#' rownames as "\code{colnames(mat)}".
#' @param Size of cells. Default autodetect. This is used to get square plot.
#' @param saveplot Save plot to given file "a.pdf", "b.png".
#' @param ... Additional parameters given to \code{\link[pheatmap]{pheatmap}} like "width", "height", .etc.
#'
#'
#' @return Nothing
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' clusterSamplePheatmap(df)
#'
clusterSamplePheatmap <- function (mat, method="pearson", digits=4,
                                   sampleAnno=NA, cellsize=NULL,
                                   cor_file=NULL, saveplot=NA, ...){
  print("Performing sample clustering using pheatmap")

  hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(100)

  pearson_cor_hc <- Matrix2colCorrelation(mat, method, digits, cor_file)
  pearson_cor <-  pearson_cor_hc$pearson_cor
  hc_result <- as.dist(1-pearson_cor)

  if(any(is.na(hc_result))){
    hc_result <- dist(t(mat))
  }

  nSample <- nrow(pearson_cor)

  if(is.null(cellsize)){
    cellsize <- 320 / nSample
  }

  pheatmap::pheatmap(pearson_cor, cluster_rows = TRUE, cluster_cols = TRUE,
                     clustering_distance_cols = hc_result,
                     clustering_distance_rows = hc_result,
                     show_colnames = T, show_rownames = T,
                     cellwidth = cellsize, cellheight=cellsize,
                     color=hmcol, annotation_row = sampleAnno,
                     annotation_col = sampleAnno, filename=saveplot)
}


# enrichMentPlot <- function(data, x, y, color_v=NULL, logGivenColumn=NULL,
#                            x_level=NULL, y_level=NULL, xval_type=NULL, sample=NULL,
#                            term_column=NULL){
#
#   #options(scipen=999)
#   # 默认sample为x轴
#   if(is.null(sample)){
#     sample = x
#   }
#
#   # 默认go term为Y轴
#   if(is.null(term_column)){
#     term_column = y
#   }
#
# xval_type = "string"
#
# if (numCheck(data[[x]])) {
# 	xval_type = "numeric"
# 	data[[x]] = mixedToFloat(data[[x]])
# }
#
# # First order by Term, then order by Sample
# if (xval_type != "numeric") {
#   data[[x]] <- factor(data[[x]], levels=x_level, ordered=T)
# 	data <- data[order(data[[term_column]], data[[sample]]), ]
# }
#
#
# if (!is.null(logGivenColumn)){
# 	log_name = paste0("negLog10_", logGivenColumn)
# 	col_name_data <- colnames(data)
# 	col_name_data <- c(col_name_data, log_name)
# 	if (! numCheck(data[[logGivenColumn]])) {
# 		stop("logGivenColumn column is <strong>not</strong> <mark>numerical</mark> column. Plase do <strong>not</strong> set log10 transform on this column.\n")
# 	} else {
# 		data[[logGivenColumn]] = mixedToFloat(data[[logGivenColumn]])
# 	}
# 	data$log_name <- log10(data[[logGivenColumn]]) * (-1)
# 	data$log_name[data$log_name==Inf] = max(data$log_name[data$log_name!=Inf]) + 2
# 	colnames(data) <- col_name_data
# 	color_v = log_name
# }
#
# if (! numCheck(data[[color_v]])) {
# 	stop("<strong>Color</strong> variable must be <mark>numbers</mark>.")
# }
#
# data[[color_v]] = mixedToFloat(data[[color_v]])
#
# # Get the count of each unique Term
# data_freq <- as.data.frame(table(data[[term_column]]))
#
# colnames(data_freq) <- c(term_column, "IDctct")
#
# data2 <- merge(data, data_freq, by=term_column)
#
# # 如果sample列是非数字
# if (!numCheck(data[[sample]])){
# 	# Collapse sample for each Term
# 	data_samp <- ddply(data2, term_column, summarize,
# 		sam_ct_ct_ct=paste(sample, collapse="_"))
#
# 	data2 <- merge(data2, data_samp, by=term_column)
#
# 	#print(data2)
#
# 	data3 <- data2[order(data2$IDctct, data2$sam_ct_ct_ct, data2[[sample]], data2[[x]],
# 	                     data2[[color_v]]), ]
# 	}
# } else{
# 	if (xval_type != "string"){
# 		data3 <- data2[order(data2$IDctct, data2$SampleGroup, data2$negLog10_qvalue), ]
# 	} else {
# 		data3 <- data2[order(data2$IDctct, data2$negLog10_qvalue), ]
# 	}
# }
# #print(data3)
#
# term_order <- unique(data3$Description)
#
# data$Description <- factor(data$Description, levels=term_order, ordered=T)
#
#
#
# #print(data)
# rm(data_freq, data2, data3)
#
# if ("" != "") {
# 	data$SampleGroup <- as.factor(data$SampleGroup)
# 	shape_level <- length(unique(data$SampleGroup))
# 	shapes = (1:shape_level)%%30
# 	shape_order <- c()
#
# 	if (length(shape_order) > 1) {
# 		data$SampleGroup <- factor(data$SampleGroup, levels=shape_order, ordered=T)
# 	} else {
# 		data$SampleGroup <- factor(data$SampleGroup)
# 	}
#
# }
#
#
#
#
#
#
# color_v <- c("green", "red")
#
# p <- ggplot(data, aes(x=SampleGroup,y=Description)) + labs(x="", y="") + labs(title="")
#
# if (("Count" != "") && ("negLog10_qvalue" != "")) {
# 	p <- p + geom_point(aes(size=Count, color=negLog10_qvalue )) + 	scale_colour_gradient(low=color_v[1], high=color_v[2], name="negLog10_qvalue")
# } else if ("Count" != "") {
# 	p <- p + geom_point(aes(size=Count ))
# } else if ("negLog10_qvalue" != "") {
# 	p <- p + geom_point(aes(color=negLog10_qvalue )) + 	scale_colour_gradient(low="color_v[1]", high=color_v[2], name="negLog10_qvalue")
# }
#
# if (("" != "") && shape_level > 6) {
# 	p <- p + scale_shape_manual(values=shapes)
# }
#
# p <- p
#
# p <- p + scale_y_discrete(labels=function(x) str_wrap(x, width=60))
#
# p <- p
#
#
# p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
# if (0 != 0){
# 	p <- p +
# 	theme(axis.text.x=element_text(angle=0,hjust=0.5,
# 	vjust=1))
# }
#
# top='top'
# bottom='bottom'
# left='left'
# right='right'
# none='none'
# legend_pos_par <- right
#
#
# uwid = 0
# vhig = 0
#
# if (uwid == 0 || vhig == 0) {
# 	x_len = length(unique(data$Description))
# 	if(x_len<10){
# 		vhig = 11
# 	} else if(x_len<20) {
# 		vhig = 11 + (x_len-10)/3
# 	} else if(x_len<100) {
# 		vhig = 14 + (x_len-20)/5
# 	} else {
# 		vhig = 40
# 	}
# 	uwid = vhig
# 	if(legend_pos_par %in% c("left", "right")){
# 		uwid = 1.5 * uwid
# 	}
# }
#
# p <- p + theme(legend.position=legend_pos_par)
#
# p <- p + theme(	panel.grid = element_blank(), panel.border=element_blank(),
# 	legend.background = element_blank(),
# 	axis.line.x=element_line(size=0.4, colour="black", linetype='solid'),
# 	axis.line.y=element_line(size=0.4, colour="black", linetype='solid'),
# 	axis.ticks = element_line(size=0.4)
# 	)
#
#
# ggsave(p, filename="/var/www/html/ImageGP/Public/source/GOenrichmentplot/1544756231.txt.scatterplot.dv.pdf", dpi=300, width=uwid,
# height=vhig, units=c("cm"))
#
# }



