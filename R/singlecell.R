
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'



#' Get stacked violin plot for seurat object
#'
#' @param object Seurat object
#' @param features A vector of genes to plot
#' @param slot_name default scale.data, accept data, count
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' stackVlnSeuratPlot(object = pbmc, features = top10$gene[c(1,3,5)])
#'
stackVlnSeuratPlot <- function(object, features, slot_name="scale.data"){
  feature_expr <- as.data.frame(t(as.data.frame(slot(object@assays$RNA, slot_name))[features,]))
  feature_expr$Cluster <- as.vector(object@active.ident)
  feature_expr <- reshape2::melt(feature_expr, id.vars=c("Cluster"),
                                 variable.name="Gene", value.name="Expr")
  p <- stackVlnPlot(feature_expr, x="Cluster", y="Expr", facets="Gene")
  p
}



#' Get stack violin plot for a normal matrix
#'
#' @param data
#'
#' At least three columns needed.
#'
#' ```
#' Gene Expr  Cluster
#' Sox2 2 1
#' Sox2 1.5 1
#' Sox2 1.2 1
#' Sox2 1.2 1
#' Sox2 20 2
#' Sox2 21 2
#' Sox2 22 2
#' Sox2 23 2
#' Sox2 0.4 3
#' Sox2 0.2 3
#' Sox3 2 1
#' Sox3 2 2
#' Sox3 2 3
#'
#' ```
#'
#' @inheritParams ggplot2::aes
#' @inheritParams ggplot2::facet_wrap
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' random_v <- c(rnorm(10, mean=1, sd=0.1), rnorm(10, mean=5), rnorm(20, mean=10),
#'               rnorm(10, mean=10), rnorm(10, mean=0.2, sd=0.01), rnorm(20, mean=1))
#' data <- data.frame(Gene=c(paste0('SOX', rep(2,40)), paste0('SOX', rep(3,40))),
#'                    Expr=random_v, Cluster=rep(c(rep(1,10), rep(2,10),rep(3,20)),2))
#' stackVlnPlot(data, x="Cluster", y="Expr", facets="Gene")
#'
stackVlnPlot <- function(data, x, y, facets, fill=NULL){
  if(is.null(fill)){
    fill = x
  }
  data[[x]] <- as.factor(data[[x]])
  p <- ggplot(data, aes_string(x=x,y=y, fill=fill)) +
    geom_violin(scale="width") +
    facet_wrap(facets, ncol=1, scales="free_y", strip.position = "left",
               labeller = as_labeller(unique(data[facets]))) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          # Customize theme so that is black & white style as requested
          panel.background = element_rect(fill = NA, colour = 'black'),
          panel.grid = element_blank()) +
    ylab("")
  p
}







