

#' Generating venDiagram plot
#'
#' @param data Data file (the first column is the name of genes or other things one want to compare, the second column is set name of each items, tab separated).
#' @param item_variable Specify the column containing all items (one of column names of data).
#' @param set_variable Specify the column containing set names (one of column names of data).
#' @param select_set_to_show Specific sets to look at.
#' @param doWeights Whether to use weights to represent data sets
#' @param type Represent data sets in different shapes. Do weighted select NO if the results are wrong
#' For 2-set Venn diagrams: circles, squares.
#' For 3-set Venn diagrams: circles, squares, ChowRuskey, triangles, AWFE.
#' For 4-set Venn diagrams: ChowRuskey, AWFE, squares or ellipses.
#' For Venn diagrams on more than four sets: classic(up to 8 sets), battle(up to 9 sets).
#' @param SetLabels Whether to plot the names of the Sets. Default TRUE.
#' @param Faces If Faces = TRUE, the sets will be filled with colors.
#' @param Sets If Sets = TRUE, the boundaries of the Sets are shown.
#' @param FaceText FaceText is a character vector which may contain any of c('weight','signature','sets','elements').
#' @inheritParams base_plot_save
#' @param ... Other parameters given to `base_plot_save`
#'
#' @return A pdf file.
#' @export
#'
#' @examples
#'
#' ## Not run:
#' vennDiagram_data = "vennDiagram.data"
#' sp_vennDiagram3(data=vennDiagram_data, header = TRUE, item_variable = "Gene",  set_variable = "Sample",
#' select_set_to_show = c("Set1","Set2","Set3"), doWeights = TRUE, type = "AWFE")
#' ## End(Not run)
#'
sp_vennDiagram3 <- function (data,
                             header = TRUE,
                             item_variable = NULL,
                             set_variable = NULL,
                             select_set_to_show = NULL,
                             doWeights = FALSE,
                             type = "ellipses",
                             SetLabels = TRUE,
                             Faces = TRUE,
                             Sets = TRUE,
                             FaceText = "weight",
                             saveplot = NULL,
                             debug = FALSE,
                             saveppt = FALSE,
                             ...) {

  # remotes::install_github("js229/Vennerable")
  options(stringsAsFactors = FALSE)

  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }

  if ("character" %in% class(data)) {
    data <- sp_readTable(data, row.names = NULL, header = header)
  } else if (!"data.frame" %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  data_colnames <- colnames(data)
  if (sp.is.null(item_variable) || sp.is.null(set_variable)) {
    stop("<item_variable> and <set_variable> should be supplied.")
  }

  if (!(item_variable %in% data_colnames &&
        set_variable %in% data_colnames)) {
    stop(paste(
      set_variable,
      'or',
      item_variable,
      'must be one of column names of data!'
    ))
  }

  data <-
    data[!is.na(data[[item_variable]]), c(item_variable, set_variable)]
  data[[set_variable]] <- make.names(data[[set_variable]])
  data <- unique(data)

  labelRecord = unique(data[[set_variable]])

  num = length(labelRecord)

  generateVennData <- function(label, data) {
    label = make.names(label)
    label_content <-
      data[data[[set_variable]] == label, item_variable]
  }

  labelContent = lapply(labelRecord, generateVennData, data)

  names(labelContent) <- labelRecord

  # 移除log
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


  if (!sp.is.null(select_set_to_show)) {
    labelContent <- labelContent[select_set_to_show]
  } else if (type %in% c("circles", "squares", "ellipses")) {
    labelContent <- labelContent[1:4]
  }
  if (length(labelContent) <= 3 && type == "ellipses") {
    type = "circles"
  }
  vennplot <- Vennerable::Venn(labelContent)

  #venncolor <- compute.Venn(vennplot)
  #gp <- VennThemes(venncolor)

  #return(gp)

  #gp[["Face"]][["11"]]$fill <-"blue"
  #gp[["Face"]][["01"]]$fill <-"green"
  #gp[["Face"]][["10"]]$fill <-"red"


  if (!is.null(saveplot)) {
    base_plot_save(saveplot, ...)
  }

  Vennerable::plot(
    vennplot,
    doWeights = doWeights,
    type = type,
    show = list(
      SetLabels = SetLabels,
      Faces = Faces,
      Sets = Sets,
      FaceText = FaceText
    )
  )
  # plot(Vstem, doWeights = FALSE)


  if (!is.null(saveplot)) {
    dev.off()
  }

  if (saveppt) {
    Vennerable::plot(
      vennplot,
      doWeights = doWeights,
      type = type,
      show = list(
        SetLabels = SetLabels,
        Faces = Faces,
        Sets = Sets,
        FaceText = FaceText
      )
    )
    eoffice::topptx(filename = paste0(saveplot, ".pptx"))
    dev.off()
  }
}
