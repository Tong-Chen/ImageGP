# www.ehbio.com/Training
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'

#' Title
#'
#' @param rlogMat xx
#' @param sample xx
#' @param ehbio_output_prefix xx
#'
#' @return xx
#' @export
#'
#' @examples
#'
#' pca_run(data)
#'
pca_run <- function(rlogMat, sample, ehbio_output_prefix){

  #print("PCA analysis")
  formulaV <- c("conditions")
  #
  topn = 5000
  rlogMat_nrow = nrow(rlogMat)
  if (topn > rlogMat_nrow){
    topn = rlogMat_nrow
  }
  #
  pca_mat = rlogMat[1:topn,]
  pca_mat <- as.data.frame(t(pca_mat))
  #
  pca <- prcomp(pca_mat, scale=T)
  #
  pca_x = pca$x
  #
  pca_individual = data.frame(samp=rownames(pca_x), pca_x, sample)
  #
  write.table(pca_individual, file=paste0(ehbio_output_prefix,".DESeq2.pca_individuals.xls"), sep="\t", quote=F, row.names=F, col.names=T)
  #
  pca_percentvar <- formatC(pca$sdev^2 * 100 / sum( pca$sdev^2))
  #
  #
  if (length(formulaV)==1) {
    p <- ggplot(pca_individual, aes(PC1, PC2, color=conditions))
  } else if (length(formulaV==2)) {
    p <- ggplot(pca_data, aes(PC1, PC2, color=conditions,
                              shape=conditions))
  }
  #
  p = p + geom_point(size=3) +
    xlab(paste0("PC1: ", pca_percentvar[1], "% variance")) +
    ylab(paste0("PC2: ", pca_percentvar[2], "% variance")) +
    geom_text_repel(aes(label=samp), show.legend=F) +
    theme_classic() +
    theme(legend.position="top", legend.title=element_blank())
  #
  p
  ggsave(p, filename=paste0(ehbio_output_prefix,".DESeq2.normalized.rlog.pca.pdf"),width=13.5,height=15,units=c("cm"))
  #
  pca_percentvar <- data.frame(PC=colnames(pca_x), Variance=pca_percentvar)
  write.table(pca_percentvar, file=paste0(ehbio_output_prefix,".DESeq2.pca_pc_weights.xls"), sep="\t", quote=F, row.names=F, col.names=T)
}
