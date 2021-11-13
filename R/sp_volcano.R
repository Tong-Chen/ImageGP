
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'

#' Generating volcano plot
#'
#' @param data Data frame or data file (with header line, the first column will not be treated as the rowname, tab separated)
#' @param geneL Data frame or data file (with header line, the first column will not be treated as the rowname, tab separated)
#' @param log2fc_var Name of fold change column.
#' @param fdr_var Name of FDR or p-value column.
#' @param coordinate_flip Default FALSE meaning `log2fc_var` in `X-axis`. Specify `TRUE` to make `fdr_var` in `X-axis`.
#' @param status_col_var Name of the column containining labels of gene expression status like `"No differnece", "UP-regulated", "Down-regulated"`. If exists, this column would be used to color each group of points. This parameter has higher priority than `significance_threshold`.
#' @param significance_threshold Set the threshold for defining DE genes in format like `c(pvalue,abs_log_fodchange)`. If this parameter is given, the program will generate `status_col_var` based on given thresholds.
#' @param status_col_var_order Changing the order of status column values.
#' Normally, the unique values of status column would be sorted alphabetically. For example, the order of `"No differnece", "UP-regulated", "Down-regulated"` would be `"Down-regulated", "No differnece", "UP-regulated"`. Here we can specify their order manually like `"UP-regulated", "Down-regulated", "No differnece"` to match the color specified below.
#' @param point_color_vector Color vector. Defgault 'red, green, grey' for 'up, dw, nodiff'
#' when `status_col_var` is not given. The order of color vector must match the order
#' of levels of `status_col_var`.
#' @param log10_transform_fdr Get `-log10(pvalue)` or `-log10(fdr)` for column given to `fdr_var`. Default FALSE, accept TRUE.
#' @param max_allowed_log10p Maximum allowed `-log10(pvalue)`. Default `Inf` meaning no limitation. Normally this should be set to 3 or 4 to make the output picture beautiful.
#' @param max_allowed_log2fc Maximum allowed `log2(foldchange)`. Default `Inf` meaning no limitation. Normally this should be set to 3 or 4 to make the output picture beautiful.
#' @param point_label_var Name of columns containing labels for points (representing genes, proteins or OTUs) to be labeled. Points with `NA`,`-` or `` value in this column will not be labeled.
#' @param point_size Point size. Default 0.8. Accept a number of name of one column.
#' @param alpha Transparency of points (0-1). 0: opaque; 1: transparent.
#' @param log2fc_symmetry Make coordiante axis symmetry to generate bettter visualization. Default TRUE.
#' @param x_label Xlab label.Default "Log2 fold change".
#' @param xintercept Default 'fc' to show fold change lines. This needs `significance_threshold`. Specify NULL to hide lines.
#' @param yintercept Default 'fdr' to show fdr lines. This needs `significance_threshold`. Specify NULL to hide lines.
#' @param y_label Ylab label.Default "Negative log10 transformed qvalue".
#' @inheritParams sp_ggplot_add_vline_hline
#' @inheritParams sp_ggplot_layout
#' @param ... Parametes given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
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
#' volcanoPlot(res_output, )
#'
#' ## Not run:
#' input <- "KO-OE_all.txt"
#' volcano_plot(input,log2fc_var='Log2FoldChange',fdr_var='Padj',status_col_var='',
#'              title="sd",label="Label",log10_transform_fdr=TRUE,point_size=5)
#' ## End(Not run)


sp_volcano_plot <-
  function(data,
           geneL=NULL,
           log2fc_var="log2FoldChange",
           fdr_var="padj",
           coordinate_flip = FALSE,
           status_col_var = NULL,
           significance_threshold = c(0.05,1),
           status_col_var_order = NULL,
           point_color_vector = c("#E6AAAA", "#9AD7EA", "#E6E5E5"),
           log10_transform_fdr = TRUE,
           max_allowed_log10p = Inf,
           max_allowed_log2fc = Inf,
           title=NULL,
           point_label_var = NULL,
           log2fc_symmetry = TRUE,
           alpha = NA,
           point_size = NA,
           extra_ggplot2_cmd = NULL,
           filename = NULL,
           xtics_angle = 0,
           x_label = 'Log2 fold change',
           y_label = 'Negative log10 transformed qvalue',
           legend.position = "top",
           xintercept = 'fc',
           yintercept = 'fdr',
           ...) {
    options(warn = -1)
    options(scipen = 999)

    if (class(data) == "character") {
      data <- sp_readTable(data, row.names = NULL, stringsAsFactors=F)
    } else if (class(data) != "data.frame") {
      stop("Unknown input format for `data` parameter.")
    }



    data_colnames <- colnames(data)

    if(! (log2fc_var %in% data_colnames && fdr_var %in% data_colnames)){
      stop(paste(log2fc_var,'or',fdr_var,'must be column names of data!'))
    }

    if (!numCheck(data[[log2fc_var]])) {
      stop("Must specify log2 transformed fold change column for Fold change column.")
    }

    if (!numCheck(data[[fdr_var]])) {
      stop("Must specify p-value, padjust or fdr column for Statistical significance column.")
    }

    if (sp.is.null(status_col_var) && length(significance_threshold) == 2) {
      status_col_var <- 'DE_genes'
      self_compute_status <- TRUE
    } else {
      self_compute_status = FALSE
    }

    if(length(significance_threshold)==2){
      fdr = significance_threshold[1]
      fc = significance_threshold[2]
    }

    if (!self_compute_status) {
      if (!sp.is.null(status_col_var_order)) {
        sig_level <- status_col_var_order
      } else{
        sig_level = c()
      }
    } else{
      sig_level <- c("UP", "DW", "NoDiff")
      data[[status_col_var]] <- ifelse(data[[fdr_var]] <= fdr,
                                   ifelse(data[[log2fc_var]] >= fc, "UP",
                                          ifelse(data[[log2fc_var]] <= fc *
                                                   (-1),
                                                 "DW", "NoDiff")) , "NoDiff")
    }

    if (length(sig_level) > 1 &&
        status_col_var != log2fc_var && status_col_var !=  fdr_var) {
      data[[status_col_var]] <-
        factor(data[[status_col_var]], levels = sig_level, ordered = T)
    }

    if (log10_transform_fdr) {
      data[[fdr_var]] <- (-1) * log10(data[[fdr_var]])
      fdr <- -1*(log10(fdr))
    }

    data[[fdr_var]] <- sapply(data[[fdr_var]], function(x) {
      if (x > max_allowed_log10p)
        max_allowed_log10p + log10(abs(x))/3
      else
        x
    })

    data[[log2fc_var]] <- sapply(data[[log2fc_var]], function(x) {
      if (x > max_allowed_log2fc)
        max_allowed_log2fc + log10(abs(x))/3
      else
        x
    })

    data[[log2fc_var]] <- sapply(data[[log2fc_var]], function(x) {
      if (x < (-1)* max_allowed_log2fc)
        (-1)*max_allowed_log2fc - log10(abs(x))/3
      else
        x
    })

    # data[[fdr_var]] <-
    #   replace(data[[fdr_var]],
    #           data[[fdr_var]] > max_allowed_log10p,
    #           max_allowed_log10p + log10(max_allowed_log10p))
    #
    # data[[log2fc_var]] <-
    #   replace(data[[log2fc_var]],
    #           data[[log2fc_var]] > max_allowed_log2fc,
    #           max_allowed_log2fc + log10(max_allowed_log2fc))
    #
    # data[[log2fc_var]] <-
    #   replace(
    #     data[[log2fc_var]],
    #     data[[log2fc_var]] < (-1) * max_allowed_log2fc,
    #     (-1) * max_allowed_log2fc - log10(max_allowed_log2fc)
    #   )


    # Below codes can be referenced if one want to extend the function
    # of labeling points.
    # This requires the gene name columns should be used as rownames or
    # at specific positions.
    # Currently not suitable for online tools. So we use the simple version.

    # if(point_label_var != "CTctctCT"){
    #   if (length(label) == 1 & is.numeric(label)) {
    #     # Label top
    #     data_line$id__ct <- rownames(data_line)
    #     data_line[(label+1):(row_num-label),"id__ct"] <- NA
    #   } else {
    #     if(is.null(names(label))){
    #       # Label given ids
    #       data_line$id__ct <- NA
    #       data_line$id__ct <- label[match(rownames(data_line), label)]
    #     } else {
    #       # Label substitute ids
    #       data_line$id__ct <- NA
    #       data_line$id__ct <- label[match(rownames(data_line), names(label))]
    #     }
    #   }
    # }

    # http://zevross.com/blog/2018/09/11/writing-efficient-and-streamlined-r-code-with-help-from-the-new-rlang-package/
    log2fc_var_en = sym(log2fc_var)
    fdr_var_en    = sym(fdr_var)
    status_col_var_en = sym(status_col_var)
    p <-
      ggplot(data = data, aes(x = !!log2fc_var_en, y = !!fdr_var_en,
                               color = !!status_col_var_en))

    if(!is.na(point_size)){
      if(is.numeric(point_size)){
        p <- p + geom_point(size=point_size)
      }else{
        p <- p + geom_point() + aes(size = !!sym(point_size))
      }
    } else {
      p <- p + geom_point()
    }

    p <- p + scale_color_manual(values = point_color_vector)


    if (!is.na(alpha)) {
      p <- p + aes(alpha = alpha) + guides(alpha = F)
    }


    if (log2fc_symmetry) {
      boundary <- ceiling(max(abs(data[, log2fc_var])))
      p <- p + xlim(-1 * boundary, boundary)
    }

    data.l <- NULL

    if(!sp.is.null(geneL)){
      if (class(geneL) == "character") {
        geneL <- sp_readTable(geneL, row.names = NULL, stringsAsFactors=F)
      } else if (class(geneL) != "data.frame") {
        stop("Unknown input format for `geneL` parameter.")
      }
      if(length(geneL) == 1){
        matched_column <-
          get_matched_columns_based_on_value(data, geneL,
                                             only_allow_one_match = T)
        #print(matched_column)
        #print(head(data[matched_column[1]]))
        #print(head(geneL[matched_column[2]]))
        #print(table(data[matched_column[1]] %in% geneL[matched_column[2]]))
        data.l <- data[data[[matched_column[1]]] %in% geneL[[matched_column[2]]],]
        point_label_var <- matched_column[1]
      } else {
        data.l <- merge_data_with_auto_matched_column(geneL, data, suffixes=c('','.y'))
        if(sp.is.null(point_label_var)){
          point_label_var <- colnames(geneL)[1]
        }
      }

      #print(head(data.l))
      #print(point_label_var)
    }

    if (!sp.is.null(point_label_var)) {
      # Only kept for back-compatible
      # giving geneL is the recommended way
      if(sp.is.null(data.l)){
        data.l <- data[data[[point_label_var]] != "-" & data[[point_label_var]] != "" & data[[point_label_var]] != "NA", ,drop=F]
      }
      #print(head(data.l))
      if (dim(data.l)[1]>0){
        #print("here")
        label_en = sym(point_label_var)
        checkAndInstallPackages(list(packages1=c("ggrepel")))
        suppressPackageStartupMessages(library(ggrepel))
        p <-
          p + geom_text_repel(
            data = data.l,
            mapping = aes(
              x = !!log2fc_var_en,
              y = !!fdr_var_en,
              label = !!label_en
            ),
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines"),
            max.overlaps = 200,
            colour = "black",
            show.legend = F
          )
      }
    }
    if(!sp.is.null(xintercept)) {
      if(xintercept=='fc'){
        if(length(significance_threshold) ==  2){
          xintercept = c(-1*fc,fc)
        }else{
          xintercept = NULL
        }
      } else {
        xintercept = as.numeric(sp_string2vector(xintercept))
      }
    }
    if(!sp.is.null(yintercept)){
      if(yintercept=='fdr'){
        if(length(significance_threshold) ==  2){
          yintercept = c(fdr)
        } else {
          yintercept = NULL
        }
      } else {
        yintercept = as.numeric(sp_string2vector(yintercept))
      }
    }

    p <- sp_ggplot_add_vline_hline(
      p,
      custom_vline_x_position = xintercept,
      custom_hline_y_position = yintercept
    )

    p <- sp_ggplot_layout(p,
                          filename = filename,
                          xtics_angle = xtics_angle,
                          legend.position = legend.position,
                          extra_ggplot2_cmd = extra_ggplot2_cmd,
                          x_label = x_label,
                          y_label = y_label,
                          title = title,
                          coordinate_flip = coordinate_flip,
                  ...)
    return(p)

  }


