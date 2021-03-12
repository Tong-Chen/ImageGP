## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ImageGP)
library(plotrix)
library(RColorBrewer)

set_data = "Set.data"

flower_plot(set_data)
# flower_plot(set_data, saveplot="Set.data.flower.pdf")

## -----------------------------------------------------------------------------
library(ImageGP)
library(ggplot2)
library(ggrepel)

res_output <- data.frame(log2FoldChange=rnorm(3000), row.names=paste0("ImageGP",1:3000))
res_output$padj <- 20 ^ (-1*(res_output$log2FoldChange^2))

padj = 0.05
log2FC = 1
res_output$level <- ifelse(res_output$padj<=padj,
                           ifelse(res_output$log2FoldChange>=log2FC,
                                  paste("groupA","UP"),
                                  ifelse(res_output$log2FoldChange<=(-1)*(log2FC),
                                         paste("groupB","UP"), "NoDiff")) , "NoDiff")
head(res_output)

# data=res_output;
# log2fc_var="log2FoldChange";
# fdr_var="padj";
#            coordinate_flip = FALSE;
#            status_col = NULL;
#            significance_threshold = c(0.05, 1);
#            status_col_level = c();
#            point_color_vector = c("red", "green", "grey");
#            log10_transform_fdr = TRUE;
#            max_allowed_log10p = Inf;
#            title='';
#            point_label_var = 'CTctctCT';
#            log2fc_symmetry = TRUE;
#            alpha = NA;
#            point_size = 0.8;
#            extra_ggplot2_cmd = NULL;
#            file_name = NULL;
#            xtics_angle = 0;
#            x_label = 'Log2 fold change';
#            y_label = 'Negative log10 transformed qvalue';
#            legend.position = "right"

sp_volcano_plot(data=res_output,
           log2fc_var="log2FoldChange",
           fdr_var="padj",
           status_col = 'level'
           )


## -----------------------------------------------------------------------------
sp_volcano_plot(data=res_output,
           log2fc_var="log2FoldChange",
           fdr_var="padj",
           significance_threshold = c(0.05, 1)
           )

## -----------------------------------------------------------------------------
# Generate one column containing genes to be labels with their symbol.
# One can also create this column easily using Excel.

label = c("Pou5f1","Sox2")
names(label) = c("ImageGP1","ImageGP4")
res_output$Symbol <- label[match(rownames(res_output), names(label))]
head(res_output)

sp_volcano_plot(data=res_output,
           log2fc_var="log2FoldChange",
           fdr_var="padj",
           status_col = 'level',
           point_label_var = "Symbol"
           )

## -----------------------------------------------------------------------------
library(ImageGP)
library(ggplot2)
library(dplyr)

manhattan_data = "manhattan.data"

sp_manhattan2_plot(data=manhattan_data, ID_var='ID', FDR_var='FDR', title="test1", point_size=2, point_label_var = "Labels")

