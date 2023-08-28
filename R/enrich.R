# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'

#' KEGG enrichment for model organism.
#'
#' @param Two columns file with first column containing DE genes and second column with group names.
#' @param org_db R annotation package like org.Hs.eg.db.
#' @param output_prefix Output prefix.
#' @param organism KEGG supported organisms like human, mouse.
#' @param pvalueCutoff Default 0.05
#' @param qvalueCutoff Default 0.2
#' @param setReadable Transfer gene ids to gene symbol. Default True.
#'
#' @return A dataframe containing top 10 enriched terms.
#' @export
#'
#' @examples
#' enrichKEGG_model(de_file, output_prefix=output_prefix)
enrichKEGG_model <-
  function(de_file,
           org_db = "org.Hs.eg.db",
           output_prefix = NULL,
           organism = "human",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           setReadable = TRUE) {
    suppressMessages(library(org_db, character.only = T))
    suppressMessages(library(DOSE))
    library(clusterProfiler)

    data <- read.table(de_file,
                       sep = "\t",
                       comment = "",
                       quote = "")
    colnames(data) <- c('gene', 'samp')
    sampC <- unique(data$samp)

    all_result <- list()

    for (samp in sampC) {
      id <- unique(data[data$samp == samp, 1])
      print(paste0("KEGG enrichment for ", samp))

      kk <- enrichKEGG(
        id,
        organism = organism,
        keyType = 'kegg',
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH",
        qvalueCutoff = qvalueCutoff
      )
      if (is.null(kk)) {
        print(paste0("No KEGG enrichemnt result for ", samp, "!!!"))
        next
      }
      if (setReadable) {
        result <- as.data.frame(setReadable(kk,
                                            org_db, keyType = "ENTREZID"))
      } else {
        result <- kk@result
      }

      if (nrow(result) > 1) {
        output <- paste(output_prefix, samp, "KEGG.xls", sep = ".")
        result$Group <- samp

        write.table(
          result,
          file = output,
          quote = F,
          sep = "\t",
          row.names = F,
          col.names = T
        )
        num = 10
        if (num > nrow(result)) {
          num = nrow(result)
        }
        all_result[[samp]] = result[1:num, ]
      }
    }

    output <- paste(output_prefix, "all.KEGG.xls", sep = ".")
    result <- do.call(rbind, all_result)
    write.table(
      result,
      file = output,
      quote = F,
      sep = "\t",
      row.names = F,
      col.names = T
    )
    invisible(all_result)
  }

#' GO enrichment for model organism.
#'
#' @param de_file Two columns file with first column containing DE genes and second column with group names.
#' @param org_db R annotation package like org.Hs.eg.db.
#' @param output_prefix Output prefix.
#' @param pvalueCutoff Default 0.05
#' @param qvalueCutoff Default 0.2
#' @param setReadable Transfer gene ids to gene symbol. Default True.
#' @param typeL A vector with default ad c("BP", "MF", "CC").
#'
#' @return NULL
#' @export
#'
#' @examples
#' enrichGO_model(de_file, output_prefix=output_prefix)
enrichGO_model <-
  function(de_file,
           org_db = "org.Hs.eg.db",
           output_prefix = NULL,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           setReadable = TRUE,
           typeL = c("BP", "MF", "CC")) {
    suppressMessages(library(org_db, character.only = T))
    suppressMessages(library(DOSE))
    suppressMessages(library(clusterProfiler))

    data <- read.table(de_file,
                       sep = "\t",
                       comment = "",
                       quote = "")
    colnames(data) <- c('gene', 'samp')
    sampC <- unique(data$samp)

    for (type in typeL) {
      all_result <- list()
      for (samp in sampC) {
        id <- unique(data[data$samp == samp, 1])

        print(paste0(type, " enrichment for ", samp))
        BP <- enrichGO(
          id,
          org_db,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = "BH",
          qvalueCutoff = qvalueCutoff,
          ont = type,
          readable = setReadable
        )
        if (is.null(BP)) {
          print(paste0("No ", type, " enrichemnt result for ", samp, "!!!"))
          next
        }

        result <-
          simplify(BP,
                   cutoff = 0.7,
                   by = "p.adjust",
                   select_fun = min)
        result <- as.data.frame(result)
        if (nrow(result) > 1) {
          result$Group <- samp

          output <- paste(output_prefix,
                          "entrez",
                          samp,
                          paste0(type, "_GO.xls"),
                          sep = ".")
          write.table(
            result,
            file = output,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T
          )

          num = 10
          if (num > nrow(result)) {
            num = nrow(result)
          }
          all_result[[samp]] = result[1:num, ]
        }
      }

      output <- paste(output_prefix, type, "all.xls", sep = ".")
      result <- do.call(rbind, all_result)
      write.table(
        result,
        file = output,
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = T
      )
      #invisible(all_result)
    }
  }

#' Enrichment of customrized pathway or annotation for any organism.
#'
#' @param de_file Two columns file with first column containing DE genes and second column with group names.
#' @param anno_file Two columns file with first column of annotations or pathways and second columns of genes within this annotation or pathway.
#' @param output_prefix Output prefix.
#' @param pvalueCutoff Default 0.05
#' @param qvalueCutoff Default 0.2
#' @param setReadable Transfer gene ids to gene symbol. Default True.
#' @param typeL A vector with default ad c("BP", "MF", "CC").
#'
#' @return NULL
#' @export
#'
#' @examples
#' enrichCustomizedPathway(de_file, output_prefix=output_prefix)
enrichCustomizedPathway <-
  function(de_file,
           anno_file,
           output_prefix = NULL,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) {
    self_anno <- read.table(anno_file,
                            header = T,
                            sep = "\t",
                            quote = "")

    colnames(self_anno) <- c("ont", "gene")

    data <- read.table(de_file,
                       sep = "\t",
                       comment = "",
                       quote = "")
    colnames(data) <- c('gene', 'samp')
    sampC <- unique(data$samp)

    all_result <- list()

    for (samp in sampC) {
      id <- unique(data[data$samp == samp, 1])
      print(paste0("Customized enrichment for ", samp))

      # self_enrich与之前enrichGO的输出结果格式一致
      self_enrich <-
        enricher(
          id,
          TERM2GENE = self_anno,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          qvalueCutoff = 0.2
        )

      result <- self_enrich@result
      if (nrow(result) > 1) {
        result$Group <- samp

        output <- paste(output_prefix, samp, "enriched.xls", sep = ".")
        result$Group <- samp
        write.table(
          result,
          file = output,
          quote = F,
          sep = "\t",
          row.names = F,
          col.names = T
        )

        num = 10
        if (num > nrow(result)) {
          num = nrow(result)
        }
        all_result[[samp]] = result[1:num, ]
      }
    }

    output <- paste(output_prefix, "all.enriched.xls", sep = ".")
    result <- do.call(rbind, all_result)
    write.table(
      result,
      file = output,
      quote = F,
      sep = "\t",
      row.names = F,
      col.names = T
    )

    invisible(result)
  }
