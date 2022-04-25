# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'


#' Transfer day to hours
#'
#' @param x A string like 1d,2d
#'
#' @return A number
#' @export
#'
#' @examples
#'
#' dh('1d')
#'
dh <- function(x) {
  if (str_extract_all(x, "[a-z]") == "h") {
    gsub("h", "", x)
  } else {
    x_time <- gsub("d", "", x)
    as.numeric(x_time) * 24
  }
}

#' Read in expression matrix and trait data (if possible).
#'
#' @param exprMat Gene expression matrix in format as "Genes x Samples".
#' The first column (gene names) must be unique among all rows and will be treated as
#' rownames. The first row (sample names) must be unique among all columns and will be
#' treated as colnames. Columns should be separted by "TAB".
#'
#' The expression data can be log transformed FPKM/TPM/CPM, \code{\link[DESeq2]{vst}} or
#' \code{\link[DESeq2]{rlog}} transformed value.
#'
#' ```
#' ID Samp1 Samp2 ... SampX
#' Gene1  1.5 2.0 ... 10
#' Gene2  1.2 4.0 ... 10
#' .
#' .
#' .
#' Gene3  2.5 2.0 ... 8
#' ```
#' @param traitData Sample attribte data with first column as sample names and other
#' columns as sample attributes. Specifically for categorical attributes, each
#' attribute one column, `0` represents not belong to while `1` represents belonging to. Or
#' one can give categorical attributes separately to "categoricalTrait".
#'
#' ```
#' ID      WT      KO      OE Height Weight Diameter
#' samp1   1       0       0       1       2       3
#' samp2   1       0       0       2       4       6
#' samp3   0       1       0       10      20      50
#' samp4   0       1       0       15      30      80
#' samp5   0       0       1       NA      9       8
#' samp6   0       0       1       4       8       7
#'
#' ```
#'
#' @param categoricalTrait Categorical attributes file with format described below.
#' The program will transferred it to 0-1 matrix like them in "traitData".
#' One can give only `traitData` or `categoricalTrait` or both (the program will
#' bind them together).
#' ```
#' ID group family
#' samp1 WT A
#' samp2 WT B
#' samp3 KO A
#' samp4 KO B
#' samp5 OE A
#' samp6 OE B
#' ```
#'
#' @inheritParams  utils::read.table
#' @param ... Other parameters given to \code{\link{read.table}}.
#' @return A lsit with three elements: `datExpr` and `traitData`, `traitColors`.
#' @export
#'
#' @examples
#'
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#'



WGCNA_readindata <-
  function(exprMat,
           traitData = NULL,
           categoricalTrait = NULL,
           sep = "\t",
           row.names = 1,
           header = T,
           quote = "",
           comment = "",
           check.names = F,
           ...) {
    cat(sp_current_time(), "Start reading in exprMat and traitData.\n")
    datExpr <-
      read.table(
        exprMat,
        sep = sep,
        row.names = row.names,
        header = header,
        quote = quote,
        comment = comment,
        check.names = check.names
      )

    sampleName = colnames(datExpr)

    continuous = 0
    category = 0
    traitColors = NULL
    result <- list(datExpr = datExpr,
                   traitData = NULL,
                   traitColors = NULL)
    if (!sp.is.null(traitData)) {
      traitData <-
        read.table(
          file = traitData,
          sep = sep,
          row.names = row.names,
          header = header,
          quote = quote,
          comment = comment,
          check.names = check.names
        )

      coln <- colnames(traitData)
      traitData = traitData[match(sampleName, rownames(traitData)),,drop=FALSE]
      #colnames(traitData) <- coln

      traitfile <- data.frame(row.names=rownames(traitData))
      for (name in coln) {
        # print (name)
        if (is.numeric(traitData[, name])) {
          numer_col <- as.data.frame(traitData[, name])
          colnames(numer_col) <- name
          traitfile[, name] <- numer_col
          # numer_col<- cbind(blank_col,numer_col)
        } else if (length(unique(traitData[, name])) == 1) {
          donothing = 2
        } else if (name == "Time") {
          time_col <- as.data.frame(apply(xc, 1, dh))
          rownames(time_col) <- row.names(traitData)
          colnames(time_col) <- name
          traitfile[, name] <- time_col
        } else if (name == "Titeration") {
          titer_col <-
            t(data.frame(str_extract_all(traitData[, name], "[0-9]+[0-9]")))
          rownames(titer_col) <- row.names(traitData)
          colnames(titer_col) <- name
          traitfile[, name] <- titer_col
          # time_titer_col<- cbind(blank_col,time_titer_col)
        } else {
          everycol <- model.matrix(as.formula(paste0("~ 0 +", name)), data = traitData)
          #everycol <- model.matrix( ~ 0 + get(name), data = traitData)
          #colnames(everycol) <-
          #  gsub("get\\(name\\)", name, colnames(everycol))
          traitfile <- cbind(traitfile, everycol)
        }
      }
      traitData <- traitfile
      continuous = 1
    }

    if (!sp.is.null(categoricalTrait)) {
      categoricalTrait <-
        read.table(
          file = categoricalTrait,
          sep = sep,
          row.names = row.names,
          header = header,
          quote = quote,
          comment = comment,
          check.names = check.names
        )

      categoricalTrait = categoricalTrait[match(sampleName, rownames(categoricalTrait)),]
      category = 1
    }

    if (continuous + category == 2) {
      traitData = cbind(traitData, categoricalTrait)
    } else if (category == 1) {
      traitData = categoricalTrait
    }

    if (continuous + category >= 1) {
      # Convert traits to a color representation:
      # white means low, red means high, grey means missing entry
      traitColors = WGCNA::numbers2colors(traitData, signed = FALSE)

      colnames(traitColors) <- colnames(traitData)
      rownames(traitColors) <- rownames(traitData)
    }
    result <-
      list(datExpr = datExpr,
           traitData = traitData,
           traitColors = traitColors)

    cat(sp_current_time(), "Finish reading in exprMat and traitData.\n")
    return(result)
  }


#' This is used to read in and check the distribution of WGCNA data in case there is
#' systematic shift of expression data (espacially for data from different detection
#' platforms).
#'

#' @param datExpr Normal gene expression matrix (gene x sample).
#'
#' @param ... Other parameters given to \code{\link{widedataframe2boxplot}}.
#'
#' @return A ggplo2 object.
#' @export
#'
#' @examples
#'
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#'
WGCNA_dataCheck <- function(datExpr, ...) {
  ##导入数据##
  # 格式如前面描述
  # 常规表达矩阵，log2转换后或
  # Deseq2的varianceStabilizingTransformation转换的数据
  # 如果有批次效应，需要事先移除，可使用removeBatchEffect
  # 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，
  # 需要quantile normalization

  cat(sp_current_time(), "Start checking data distribution of input matrix.\n")
  ellipse_par <- list(...)
  p = widedataframe2boxplot(datExpr, ...)
  if (!"saveplot" %in% names(ellipse_par)) {
    print(p)
  }
  cat(sp_current_time(), "End checking data distribution of input matrix.\n")
  return(p)
}


#' Filter low variance genes by given minimal \code{\link{mad}} value or keep top
#' number/percent genes with bigger variances.
#'
#' @inheritParams WGCNA_dataCheck
#' @param minimal_mad Minimal allowed mad value.
#' @param top_mad_n An integer larger than 1 will be used to get top x genes (like top 5000).
#' A float number less than 1 will be used to get top x fraction genes (like top 0.7 of
#' all genes).
#' @param rmVarZero Default TRUE. Remove genes with variance as 0. Normally for PCA or
#' correlation analysis.
#' @param noLessThan Specify the lowest number of genes to be kept. Default `NULL` meaning no lower limit.
#' @param value_type Specify the way for statistical computation. Default mad, accept mean, var.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3)
#' dataFilter(df)
#'
dataFilter <-
  function(datExpr,
           minimal_mad = NULL,
           top_mad_n = 0.75,
           rmVarZero = T,
           noLessThan = NULL,
           value_type = mad) {
    if(mode(value_type) == "character"){
      value_type = switch(
        value_type,
        mean = mean,
        var = var,
        mad = mad,
        median = median,
        sum = sum
      )
    }

    m.mad <- apply(datExpr, 1, value_type)

    if (!sp.is.null(minimal_mad)) {
      datExpr <- datExpr[which(m.mad > minimal_mad), ]
    } else {
      datExpr <- datExpr[order(m.mad, decreasing = T), ]
      nGenes <- nrow(datExpr)
      if (top_mad_n < 1) {
        top_mad_n = ceiling(nGenes * top_mad_n)
      } else if (top_mad_n == 1) {
        top_mad_n = nGenes
      }

      if (!sp.is.null(noLessThan)) {
        if (top_mad_n < noLessThan) {
          top_mad_n = noLessThan
        }
      }

      if (top_mad_n > nGenes) {
        top_mad_n = nGenes
      }
      datExpr <- datExpr[1:top_mad_n, ]
    }

    if (rmVarZero) {
      m.var <- apply(datExpr, 1, var)
      datExpr <- datExpr[which(m.var > 0), ]
    }

    return(datExpr)
  }

#' Filter low variance genes by given minimal \code{\link{mad}} or \code{\link{var}} or
#' other statistical value or keep top
#' number/percent genes with bigger variances.
#'
#' @inheritParams WGCNA_dataCheck
#' @param statistical_value_type Specify the way for statistical computation. Default mad, accept mean, var, sum, median.
#' @param minimal_threshold Minimal allowed statistical value.
#' @param top_n An integer larger than 1 will be used to get top x genes (like top 5000).
#' A float number less than 1 will be used to get top x fraction genes (like top 0.7 of
#' all genes).
#' @param rmVarZero Default TRUE. Remove genes with variance as 0. Normally for PCA or
#' correlation analysis.
#' @param noLessThan Specify the lowest number of genes to be kept. Default `NULL` meaning no lower limit.
#' @param keep_filtered_as_others Get sums of all filtered items as an new item - Others. Default FALSE.

#' @return A dataframe.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3)
#' dataFilter2(df)
#'
dataFilter2 <-
  function(datExpr,
           minimal_threshold = NULL,
           top_n = 1,
           rmVarZero = F,
           noLessThan = NULL,
           statistical_value_type = mad,
           keep_filtered_as_others = F) {
    if(top_n == 1) {
      return(datExpr)
    }

    nGenes <- nrow(datExpr)

    if (top_n >= nGenes) {
      return(datExpr)
    }

    if(mode(statistical_value_type) == "character"){
      statistical_value_type = switch(
        statistical_value_type,
        mean = mean,
        var = var,
        mad = mad,
        median = median,
        sum = sum
      )
    }

    m.mad <- apply(datExpr, 1, statistical_value_type)

    if (sp.is.null(minimal_threshold)) {
      m.mad.sorted <- sort(m.mad, decreasing = T)

      if (top_n < 1) {
        top_n = ceiling(nGenes * top_n)
      }

      if (!sp.is.null(noLessThan)) {
        if (top_n < noLessThan) {
          top_n = noLessThan
        }
      }
      print(top_n)
      minimal_threshold <- m.mad.sorted[top_n]
    }

    datExpr2 <- datExpr[which(m.mad >= minimal_threshold), ]

    if(nrow(datExpr2)==0){
      stop("ALl data are filtered. Please relax the filtering threshold.")
    }

    if(keep_filtered_as_others){
      datExprFiltered = datExpr[which(m.mad < minimal_threshold), ]
      if(nrow(datExprFiltered)>0){
        others <- colSums(datExprFiltered)
        datExpr2 <- rbind(datExpr2, others)
        rownames(datExpr2)[nrow(datExpr2)] = "Others"
      }
    }


    if (rmVarZero) {
      m.var <- apply(datExpr2, 1, var)
      datExpr2 <- datExpr2[which(m.var > 0), ]
    }

    return(datExpr2)
  }



#' Filter data for WGCNA input to increase computing efficiency without loosing too
#' many information.
#'
#' @param wgcnaL A matrix or an object return by \code{WGCNA_readindata}.
#' @param ... Parameters given to \code{dataFilter}.
#'
#' @return A dataframe (samples x genes).
#' @export
#'
#' @examples
#'
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#'
WGCNA_dataFilter <- function (wgcnaL, ...) {
  cat(sp_current_time(), "Start filtering data.\n")
  if(class(wgcnaL) == "list"){
    datExpr = wgcnaL$datExpr
  } else {
    datExpr = wgcnaL
  }
  datExpr <- dataFilter(datExpr, noLessThan = 1000, ...)

  ## 转换为样品在行，基因在列的矩阵
  datExpr <- as.data.frame(t(datExpr))

  ## 检测缺失值
  gsg = WGCNA::goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes) > 0)
      cat("\tRemoving genes:",
                  paste(names(datExpr)[!gsg$goodGenes], collapse = ","),"\n")

    if (sum(!gsg$goodSamples) > 0)
      cat("\tRemoving samples:",
                  paste(rownames(datExpr)[!gsg$goodSamples], collapse = ","),"\n")

    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }


  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

  cat(sp_current_time(), "After filtering,",
      nGenes,
      'genes x',
      nSamples,
      "samples remained.\n")
  if(class(wgcnaL) == "list"){
    wgcnaL$traitData = wgcnaL$traitData[rownames(datExpr),,drop=F]
    wgcnaL$traitColor = wgcnaL$traitColor[rownames(datExpr),,drop=F]
    wgcnaL$datExpr = datExpr
  } else {
    wgcnaL = datExpr
  }
  return(wgcnaL)
}

#' Keep samples in trait data the same of exprMat.
#'
#' @param datExpr expression matrix
#' @param trait trait data matrix
#'
#' @return balanced trait data frame
#' @export
#'
#' @examples
#'
#' traitData = WGCNA_filterTrait(datExpr, traitData)
#' traitColor = WGCNA_filterTrait(datExpr, traitColor)
#'
WGCNA_filterTrait <- function(datExpr, trait){
  return(trait[rownames(datExpr),,drop=F])
}

#' Sample cluster and outlier detection
#'
#' @param wgcnaL A matrix or an object return by \code{WGCNA_readindata}. A transformed gene expression matrix normally output by \code{WGCNA_dataFilter}.
#' Samples x Genes.
#' @param thresholdZ.k Threshold for defining outliers. First compute the overall
#' corelation of one sample to other samples. Then do Z-score transfer for all
#' correlation values. The samples with corelation values less than given value
#' would be treated as outliers.
#' Default -2.5 meaning -2.5 std.
#' @param traitColors Sample attributes data frame transferred by
#' \code{\link[WGCNA]{numbers2colors}} or generated in \code{\link{WGCNA_readindata}}.
#' @inheritParams base_plot_save
#' @param removeOutlier Remove outlier samples. Normally this should be only performed if
#' no suitable soft power can be found.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#'
#'
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#'
WGCNA_sampleClusterDetectOutlier <-
  function(wgcnaL,
           thresholdZ.k = -2.5,
           saveplot = NULL,
           removeOutlier = F,
           traitColors = NULL,
           ...) {
    cat(sp_current_time(), "Detect outlier samples.\n")
    if(class(wgcnaL) == "list"){
      datExpr = wgcnaL$datExpr
    } else {
      datExpr = wgcnaL
    }
    ## 样本层级聚类，查看有无离群值
    # sample network based on squared Euclidean distance note that we
    # transpose the data
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

    cat("\tThreshold for detecting outlier samples are", thresholdZ.k,"\n")
    # the color vector indicates outlyingness (red)
    outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")

    # calculate the cluster tree using flahsClust or hclust
    sampleTree = hclust(as.dist(1 - A), method = "average")
    # Convert traits to a color representation: where red indicates high
    # values
    if (!is.null(traitColors)) {
      # traitColors = data.frame(WGCNA::numbers2colors(datTraits, signed = FALSE))
      #dimnames(traitColors)[[2]] = paste(names(datTraits), "C", sep = "")
      datColors = data.frame(outlierC = outlierColor, traitColors)
    } else {
      datColors = data.frame(outlierC = outlierColor)
    }
    # Plot the sample dendrogram and the colors underneath.

    if (!sp.is.null(saveplot)) {
      base_plot_save(saveplot, ...)
    }

    WGCNA::plotDendroAndColors(
      sampleTree,
      groupLabels = names(datColors),
      colors = datColors,
      main = "Sample dendrogram with/without trait heatmap"
    )

    if (!sp.is.null(saveplot)) {
      dev.off()
    }

    if (removeOutlier) {
      cat("\tRemoving outlier samples <",
          rownames(datExpr[which(Z.k < thresholdZ.k), ]),">.\n")
      datExpr <- datExpr[which(Z.k >= thresholdZ.k), ]

      cat("\tAfter removing outlier samples, ", nrow(datExpr), "samples kept.\n")
    }
    if(class(wgcnaL) == "list"){
      wgcnaL$traitData = wgcnaL$traitData[rownames(datExpr),,drop=F]
      wgcnaL$traitColor = wgcnaL$traitColor[rownames(datExpr),,drop=F]
      wgcnaL$datExpr = datExpr
    } else {
      wgcnaL = datExpr
    }

    cat(sp_current_time(), "Finish detecting outlier samples.\n")
    return(wgcnaL)
  }



#' Select soft power, the minimum number to get scale Free Topology Model Fit value (R2) larger than 0.85.
#'
#' @param networkType Default "signed". Allowed values are (unique abbreviations of)
#' "unsigned", "signed", "signed hybrid". Correlation and distance are transformed as
#' follows:
#'
#' 1. for type = "unsigned", adjacency = |cor|^power;
#'
#' 2. for type = "signed", adjacency = (0.5 * (1+cor) )^power;
#'
#' 3. for type = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise;
#'
#' and for type = "distance", adjacency = (1-(dist/max(dist))^2)^power.
#'
#' @inheritParams  WGCNA_sampleClusterDetectOutlier
#' @param maxPower Specify maximum power to check. Default 30 for "unsigned" network
#' and 40 for other type. Any number less than 20 would be treated as 20.
#' @param RsquaredCut R2 for defining scale-free network (default 0.85). Any number larger than 1 would be treated as 0.99.
#'
#' @return A list
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' #datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#'
WGCNA_softpower <-
  function(datExpr,
           networkType = "signed",
           saveplot = NULL,
           maxPower = NULL,
           RsquaredCut = 0.85,
           ...) {
    ## 软阈值筛选
    ## 软阈值的筛选原则是使构建的网络更符合无标度网络特征。
    cat(sp_current_time(), "Select soft power for network construction.\n")
    if (sp.is.null(maxPower)) {
      if (networkType == "unsigned") {
        maxPower = 30
      } else {
        maxPower = 45
      }
    }
    if (maxPower < 20) {
      maxPower = 20
    }
    if (RsquaredCut >= 1) {
      RsquaredCut = 0.99
    }

    cat("\tMax power used for selection is", maxPower,".\n")

    powers = c(c(1:10), seq(from = 12, to = maxPower, by = 2))

    sft = WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = networkType,
      verbose = 5,
      RsquaredCut = RsquaredCut
    )


    power = sft$powerEstimate

    power_data = sft$fitIndices
    power_data$Choose = "Other tested powers"



    # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
    # 网络越符合无标度特征 (non-scale)



    ## ---- echo=T-------------------------------------------------------------
    # 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
    # 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
    # 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
    # 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
    # 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。

    if (is.na(power)) {
      cat("\tUsing experience power since no suitable power found.\n")

      nSamples = nrow(datExpr)
      power = ifelse(
        nSamples < 20,
        ifelse(networkType == "unsigned", 9, 18),
        ifelse(
          nSamples < 30,
          ifelse(networkType == "unsigned", 8, 16),
          ifelse(
            nSamples < 40,
            ifelse(networkType == "unsigned", 7, 14),
            ifelse(networkType == "unsigned", 6, 12)
          )
        )
      )
      power_choose = which(power_data$Power == power)
      power_data$Choose[power_choose] = "Classicial_power"
      power_text = "No suitable power found. Use experience power: "
    } else {
      power_choose = which(power_data$Power == power)
      power_data$Choose[power_choose] = "Choosed_power"
      power_text = "Finally choosed power is: "
    }

    cat("\tFinally chooosed power is :", power, ".\n")

    # 图展示在文档中

    p1 <- ggplot(power_data, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
      xlab("Soft Threshold (power)") +
      ylab("Scale Free Topology Model Fit, signed R^2") +
      ggtitle(paste0(power_text, power, ".\n", "Scale independence")) +
      geom_text(aes(label = Power, color = Choose))

    p1 <- sp_ggplot_add_vline_hline(
      p1,
      custom_hline_y_position = RsquaredCut,
      custom_hline_anno = paste0("R squared Cut ", RsquaredCut),
      color = "red"
    ) + theme_classic()
    p1 <- sp_ggplot_layout(p1, legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5))


    p2 <- ggplot(power_data, aes(x = Power, y = mean.k.)) +
      xlab("Soft Threshold (power)") +
      ylab("Mean Connectivity") +
      ggtitle("Mean connectivity") +
      geom_point(aes(color = Choose))
    p2 <-
      sp_ggplot_add_vline_hline(
        p2,
        custom_hline_y_position = power_data$mean.k.[power_choose],
        custom_hline_anno = round(power_data$mean.k.[power_choose]),
        custom_hline_anno_x_pos = power,
        color = "red"
      ) + theme_classic()
    p2 <- sp_ggplot_layout(p2, legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5))

    p = p1 %>% aplot::insert_bottom(p2)

    if (sp.is.null(saveplot)) {
      p
    } else {
      ggsave(
        p,
        filename = saveplot,
        units = c("cm"),
        width = 15,
        height = 18
      )
    }

    cat(sp_current_time(), "Finished selecting soft power for network construction.\n")
    return(list(power=power,p=p))
  }



#' WGCNA main function.
#'
#' @inheritParams WGCNA::blockwiseModules
#' @inheritParams base_plot_save
#' @param width Width of graphics in inches. Default 14.
#' @param height Height of graphics in inches. Default 7.
#' @param ... Other parameters given to \code{\link[WGCNA]{blockwiseModules}}.
#' @param noplot Return construccted net object only.
#' @param dynamicCutPlot Plot merged modules as well as dynamic cutted modules before merge.
#'
#' @return A network
#' @export
#'
#' @section Explanations:
#'
#' * maxBlockSize: Maximum allowed number of genes in one block.
#' Default "all genes in one block". Please check
#' \url{http://www.peterlangfelder.com/blockwise-network-analysis-of-large-data/} for
#' more explanations. Analyzing a set of 20,000 genes requires between 8 and 16 GB of memory;
#' 40,000 genes would increase the requirement to 32-64 GB.
#' A full network analysis of 500,000 genes would theoretically require some 7 TB of memory.
#'
#' Quote: I emphasize that the blockwise analysis creates an approximation to the
#' network that would result from a single block analysis. The approximation
#' is often very good but the modules are not quite the same.
#' If possible, I recommend running the analysis in a single block;
#' if not, use the largest blocks your computer can handle.
#'
#' * mergeCutHeight The larger the less number of merged modules. Normally 0.15-0.3.
#'
#' * networkType: "signed" is recommended. However, number of genes in modules would also
#' be less for "signed" netwrok. Check \code{WGCNA_softpower} for the detail.
#'
#' * corType: biweight mid-correlation (bicor) recommended.
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#'
WGCNA_coexprNetwork <- function(datExpr,
                                power,
                                maxBlockSize = NULL,
                                minModuleSize = NULL,
                                networkType = "signed",
                                mergeCutHeight = 0.2,
                                numericLabels = TRUE,
                                pamRespectsDendro = FALSE,
                                saveTOMs = TRUE,
                                corType = "bicor",
                                maxPOutliers = NULL,
                                loadTOM = TRUE,
                                TOMDenom = "min",
                                deepSplit = 1,
                                stabilityCriterion = "Individual fraction",
                                saveTOMFileBase = "blockwiseTOM",
                                verbose = 3,
                                randomSeed = 1117,
                                saveplot = NULL,
                                width = 14,
                                height = 7,
                                noplot = FALSE,
                                dynamicCutPlot = TRUE,
                                ...) {
  cat(sp_current_time(), "Start constructing WGCNA network.\n")
  if (sp.is.null(maxBlockSize)) {
    maxBlockSize = ncol(datExpr)
  }

  if(sp.is.null(minModuleSize)){
    minModuleSize = min(20, ncol(datExpr)/2 )
  }

  # corFnc = ifelse(corType=="pearson", cor, bicor)
  # 对二元变量，如样本性状信息计算相关性时，
  # 或基因表达严重依赖于疾病状态时，需设置下面参数
  if (sp.is.null(maxPOutliers)) {
    maxPOutliers = ifelse(corType == "pearson", 1, 0.05)
  }

  net = WGCNA::blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = maxBlockSize,
    TOMType = networkType,
    minModuleSize = minModuleSize,
    networkType = networkType,
    mergeCutHeight = mergeCutHeight,
    numericLabels = numericLabels,
    pamRespectsDendro = pamRespectsDendro,
    saveTOMs = saveTOMs,
    corType = corType,
    maxPOutliers = maxPOutliers,
    loadTOM = loadTOM,
    TOMDenom = TOMDenom,
    deepSplit = deepSplit,
    stabilityCriterion = stabilityCriterion,
    saveTOMFileBase = saveTOMFileBase,
    verbose = verbose,
    randomSeed = randomSeed
  )

  cat("\tCore parameters for constructing WGCNA network:", "\n\t\tpower =", power,
      "\n\t\tminModuleSize =", minModuleSize, "\n\t\tnetworkType =", networkType,
      "\n\t\tmergeCutHeight =", mergeCutHeight, "\n\t\tdeepSplit =", deepSplit,"\n")
  if (noplot) {
    return(net)
  }
  # 根据模块中基因数目的多少，降序排列，依次编号为 1-最大模块数。
  # **0 (grey)**表示**未**分入任何模块的基因。
  # table(net$colors)

  ## ---- echo=T, fig.cap="层级聚类树展示各个模块"---------------------------
  ## 灰色的为**未分类**到模块的基因。
  # Convert labels to colors for plotting
  moduleLabels = net$colors
  moduleColors = WGCNA::labels2colors(moduleLabels)

  # Plot the dendrogram and the module colors underneath
  # 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间

  if (!sp.is.null(saveplot)) {
    base_plot_save(saveplot,
                   width = width,
                   height = height,
                   bg = "white",
                   ...)
  }


  if (dynamicCutPlot) {
    dynamicColors <- WGCNA::labels2colors(net$unmergedColors)
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]],
      cbind(dynamicColors, moduleColors),
      c("Dynamic Tree Cut", "Module colors"),
      dendroLabels = FALSE,
      hang = 0.5,
      addGuide = TRUE,
      guideHang = 0.05
    )
  } else {
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]],
      moduleColors,
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.5,
      addGuide = TRUE,
      guideHang = 0.05
    )
  }

  if (!sp.is.null(saveplot)) {
    dev.off()
  }
  cat(sp_current_time(), "Finished constructing WGCNA network.\n")
  return(net)
}


#' Save gene-module relationships, MEs and plot module correlations.
#'
#' @param net \code{\link{WGCNA_coexprNetwork}} or \code{\link[WGCNA]{blockwiseModules}} returned WGCNA object.
#' @inheritParams WGCNA_coexprNetwork
#' @inheritParams base_plot_save
#' @param prefix prefix for output files.
#' @param ... Additional parameters given to plot output (\code{\link{pdf}}, \code{\link{png}},...) like "width", "height", .etc.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#'
#'
WGCNA_saveModuleAndMe <-
  function(net,
           datExpr,
           prefix = "ehbio",
           saveplot = NULL,
           ...) {
    ## 共表达网络结果输出
    cat(sp_current_time(), "Save WGCNA modules.\n")
    #1. 输出基因及其所在模块信息，方便对模块进行富集分析。
    #2. 输出模块的主成分信息 (ME)，代表模块整体基因表达量

    moduleLabels = net$colors
    moduleColors = WGCNA::labels2colors(moduleLabels)

    ### 基因和所在模块信息
    gene_module <-
      data.frame(ID = colnames(datExpr), module = moduleColors)
    gene_module = gene_module[order(gene_module$module), ]
    write.table(
      gene_module,
      file = paste0(prefix, ".gene_module.xls"),
      sep = "\t",
      quote = F,
      row.names = F
    )

    # module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
    MEs = net$MEs

    ### 不需要重新计算，改下列名字就好
    ### 官方教程是重新计算的，其实可以不用这么麻烦
    MEs_col = MEs
    colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(
      stringr::str_replace_all(colnames(MEs), "ME", "")
    )))
    MEs_col = orderMEs(MEs_col)

    ## 保存模块代表性信息
    MEs_colt = as.data.frame(t(MEs_col))
    colnames(MEs_colt) = rownames(datExpr)
    write.table(
      data.frame(Module=rownames(MEs_colt),MEs_colt),
      file = paste0(prefix, ".module_eipgengene.xls"),
      sep = "\t",
      quote = F,
      row.names = F
    )

    if (ncol(MEs_col) < 4) {
      return(MEs_col)
    }

    if (!sp.is.null(saveplot)) {
      base_plot_save(saveplot, bg = "white", ...)
    }
    # marDendro/marHeatmap 设置下、左、上、右的边距
    WGCNA::plotEigengeneNetworks(
      MEs_col,
      "Eigengene adjacency heatmap",
      marDendro = c(3, 3, 2, 4),
      marHeatmap = c(3, 4, 2, 2),
      plotDendrograms = T,
      xLabelsAngle = 90
    )

    if (!sp.is.null(saveplot)) {
      dev.off()
    }
    cat(sp_current_time(), "Finished saving WGCNA modules.\n")
    return(MEs_col)
  }


#' Heatmap showing correlation among MEs and traits.
#'
#' @param MEs_col Module epigenes generated in \code{\link{WGCNA_saveModuleAndMe}}.
#' @param traitData Sample attributes data frame.
#' Or the "traitData" generated in \code{\link{WGCNA_readindata}}.
#' @inheritParams WGCNA_sampleClusterDetectOutlier
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#' WGCNA_MEs_traitCorrelationHeatmap(MEs_col, traitData=wgcnaL$traitData)
#'
WGCNA_MEs_traitCorrelationHeatmap <-
  function(MEs_col, traitData, saveplot = NULL, ...) {
    cat(sp_current_time(), "Plot module trait correlation heatmap.\n")
    # if (!sp.is.null(saveplot)) {
    #   base_plot_save(saveplot, bg = "white", ...)
    # }
    if (sp.is.null(saveplot)) {
      saveplot = NA
    }
    if (sp.is.null(traitData)) {
      # MEs_colpheno = orderMEs(MEs_col)
      MEs_colpheno = MEs_col
    } else {
      traitData <- WGCNA_filterTrait(MEs_col, traitData)
      MEs_colpheno = cbind(MEs_col, traitData)
      # MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
      #ggcorr_MEs_colpheno  <- ggcorr(MEs_colpheno, hjust = 0.75, size = 5, color = "grey50", layout.exp = 1)
      #ggsave(paste0(prefix,"_ggcorr_MEs_colpheno.pdf"))
    }

    MEs_colpheno_cor <- cor(MEs_colpheno)
    pheatmap::pheatmap(MEs_colpheno_cor, filename=saveplot)

    # WGCNA::plotEigengeneNetworks(
    #   MEs_colpheno,
    #   "Eigengene adjacency heatmap",
    #   marDendro = c(3, 3, 2, 4),
    #   marHeatmap = c(3, 4, 2, 2),
    #   plotDendrograms = T,
    #   xLabelsAngle = 90
    # )
    # if (!sp.is.null(saveplot)) {
    #   dev.off()
    # }
    cat(sp_current_time(), "Finish plotting module trait correlation heatmap.\n")
  }

#' Export WGCNA resullt for cytoscape input edges and nodes.
#'
#' @inheritParams WGCNA_saveModuleAndMe
#' @inheritParams WGCNA_coexprNetwork
#' @param TOM_plot Get TOM plot and save to file given here like 'tomplot.pdf'.
#' @param fulledge Output all edges (very large). Default FALSE. Only output edges whithin modules.
#'
#' @return A list with edgeData and nodeData as two elements.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#' WGCNA_MEs_traitCorrelationHeatmap(MEs_col, traitData=wgcnaL$traitData)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#'
WGCNA_cytoscape <-
  function(net,
           power,
           datExpr,
           TOM_plot = NULL,
           prefix = "ehbio",
           fulledge = F) {
    ## 获取TOM矩阵，导出Cytoscape可用的数据方便网络图绘制
    # 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
    # 否则需要再计算一遍，比较耗费时间
    # TOM = TOMsimilarityFromExpr(datExpr, power=power, corType=corType, networkType=type)
    cat(sp_current_time(), "Start generating input files for cytoscape.\n")
    cat("\tLoad TOM similarity matrix.\n")
    load(net$TOMFiles[1], verbose = T)
    file.remove(net$TOMFiles[1])

    TOM <- as.matrix(TOM)

    dissTOM = 1 - TOM
    # Transform dissTOM with a power to make moderately strong
    # connections more visible in the heatmap
    plotTOM = dissTOM ^ power
    # Set diagonal to NA for a nicer plot
    diag(plotTOM) = NA

    moduleLabels = net$colors
    moduleColors = WGCNA::labels2colors(moduleLabels)

    # Call the plot function

    if (!sp.is.null(TOM_plot)) {
      base_plot_save(TOM_plot, width = 20, height = 20)
      TOMplot(plotTOM, net$dendrograms, moduleColors,
              main = "Network heatmap plot, all genes")
      dev.off()
    }


    probes = colnames(datExpr)
    dimnames(TOM) <- list(probes, probes)

    # Export the network into edge and node list files Cytoscape can read
    # threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
    # cytoscape中再调整
    cyt = WGCNA::exportNetworkToCytoscape(
      TOM,
      weighted = TRUE,
      threshold = 0.01,
      nodeNames = probes,
      nodeAttr = moduleColors
    )

    edgeData <- cyt$edgeData[, c(1:3)]
    colnames(edgeData) <- c("Source", "Target", "Correlation")

    if (fulledge) {
      write.table(
        edgeData,
        paste0(prefix, ".cytoscape_full_edges.txt"),
        quote = F,
        sep = "\t",
        row.names = F
      )
    }
    nodeData <- cyt$nodeData[, c(1, 3)]
    colnames(nodeData) <- c("nodeName", "Module")
    write.table(
      nodeData,
      paste0(prefix, ".cytoscape_full_nodes.txt"),
      quote = F,
      sep = "\t",
      row.names = F
    )

    Module1 <- nodeData[match(edgeData$Source, nodeData$nodeName), 2]
    Module2 <- nodeData[match(edgeData$Target, nodeData$nodeName), 2]
    # Only keep interactions within module
    edgeData <- edgeData[Module1 == Module2, ]

    write.table(
      edgeData,
      paste0(prefix, ".cytoscape_edges_within_modules.txt"),
      quote = F,
      sep = "\t",
      row.names = F
    )

    #sapply(moduleColors, WGCNA_cytoscape_each_module, edgeData=edgeData,
    #       nodeData=nodeData, prefix=prefix)

    cyt = list(edgeData = edgeData, nodeData = nodeData)
    cat(sp_current_time(), "Finish generating input files for cytoscape.\n")
    return(cyt)
  }




#' Get top x hub genes for each module.
#'
#' @param cyt A list containing two elements (edgeData and nodeData) generated by
#' \code{WGCNA_cytoscape} (specifically onle whithin module
#' interactions are kept in edgeData).
#' @param top_hub_n A number to get top x hub genes.
#' @inheritParams WGCNA_cytoscape
#'
#' @return A dataframe containing selected hub genes.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#' WGCNA_MEs_traitCorrelationHeatmap(MEs_col, traitData=wgcnaL$traitData)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#'
WGCNA_hubgene <- function(cyt,
                          top_hub_n = 20,
                          prefix = 'ehbio') {
  #library(dplyr)
  cat(sp_current_time(), "Get top hub genes for each module.\n")
  edgeData1 <- cyt$edgeData[, c(1, 2, 3)]
  edgeData2 <- cyt$edgeData[, c(2, 1, 3)]
  nodeattrib <- cyt$nodeData
  colnames(edgeData1) <- c("Node1", "Node2", "Weight")
  colnames(edgeData2) <- c("Node1", "Node2", "Weight")
  edgeData <- rbind(edgeData1, edgeData2)
  edgeData$Module1 <-
    nodeattrib[match(edgeData$Node1, nodeattrib$nodeName), 2]
  edgeData <- edgeData[edgeData$Module1!="grey",,drop=F]
  #edgeData$Module2 <- nodeattrib[match(edgeData$Node2, nodeattrib$nodeName), 2]
  #edgeData <- edgeData[edgeData$Module1==edgeData$Module2,c(1,3,4)]
  #head(edgeData)
  edgeData <- edgeData[, c(1, 3, 4)]

  nodeTotalWeight <-
    edgeData %>% dplyr::group_by(Node1, Module1) %>% dplyr::summarise(weight =
                                                                        sum(Weight))

  nodeTotalWeight <-
    nodeTotalWeight[with(nodeTotalWeight, order(Module1,-weight)), ]
  #head(nodeTotalWeight)

  nodeTotalWeightTop20 = nodeTotalWeight %>% dplyr::group_by(Module1) %>% dplyr::top_n(top_hub_n, weight)

  colnames(edgeData1) <- c("Source", "Target", "Correlation")
  hub_edgeData = edgeData1[(edgeData1$Source %in% nodeTotalWeightTop20$Node1) &
                             (edgeData1$Target %in% nodeTotalWeightTop20$Node1), ]

  write.table(
    nodeTotalWeightTop20,
    paste(prefix, "hubgenes.txt", sep = "."),
    quote = F,
    sep = "\t",
    row.names = F
  )

  write.table(
    hub_edgeData,
    paste(prefix, "hubgenes.edges.txt", sep = "."),
    quote = F,
    sep = "\t",
    row.names = F
  )

  cat(sp_current_time(), "Finish outputing top hub genes for each module.\n")
  return(nodeTotalWeightTop20)
}



#' Module-trait heatmap
#'
#' @inheritParams WGCNA_MEs_traitCorrelationHeatmap
#' @inheritParams WGCNA_coexprNetwork
#' @inheritParams WGCNA_saveModuleAndMe
#' @param angle_x Rotation angle for x-axis labels
#' @param up_color Vector of colours to use for upper triangles (which representing pearson correlations values).
#' @param down_color Vector of colours to use for lower triangles (which representing significance p-values).
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' traitData <- wgcnaL$traitData
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#' WGCNA_MEs_traitCorrelationHeatmap(MEs_col, traitData=traitData)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#' WGCNA_moduleTraitPlot(MEs_col, traitData=traitData)
#'
#'
WGCNA_moduleTraitPlot <-
  function(MEs_col,
           traitData=NULL,
           corType = "bicor",
           saveplot = NULL,
           prefix = 'ehbio',
           angle_x = 90,
           up_color = c("red", "white", "blue"),
           down_color = c("green", "white"),
           ...) {
    ### 模块与表型数据关联
    cat(sp_current_time(), "Plot correlation heatmap among all traits and modules.\n")
    if (sp.is.null(traitData)) {
      cat(sp_current_time(), "No plot since no trait data available.\n")
      return(1)
    }else{
      #print(traitData)
      traitData <- WGCNA_filterTrait(MEs_col, traitData)
      #print(trairData)
    }
    nSamples <- nrow(traitData)
    robustY = ifelse(corType == "pearson", T, F)
    if (corType == "pearson") {
      modTraitCor = cor(MEs_col, traitData, use = "p")
      modTraitP = WGCNA::corPvalueStudent(modTraitCor, nSamples)
    } else {
      modTraitCorP = WGCNA::bicorAndPvalue(MEs_col, traitData, robustY = robustY)
      modTraitCor = modTraitCorP$bicor
      modTraitP   = modTraitCorP$p
    }
    # signif表示保留几位小数
    textMatrix = paste(signif(modTraitCor, 3), "\n(", signif(modTraitP, 2), ")", sep = "")
    dim(textMatrix) = dim(modTraitCor)

    if (!sp.is.null(saveplot)) {
      base_plot_save(saveplot, bg = "white", ...)
    }
    #labeledHeatmap(
    #  Matrix = modTraitCor,
    #  xLabels = colnames(traitData),
    #  yLabels = colnames(MEs_col),
    #  cex.lab = 0.5,
    #  ySymbols = colnames(MEs_col),
    #  colorLabels = FALSE,
    #  colors = blueWhiteRed(50),
    #  textMatrix = textMatrix,
    #  setStdMargins = FALSE,
    #  cex.text = 0.5,
    #  zlim = c(-1, 1),
    #  main = paste("Module-trait relationships")
    #)
    module_name = colnames(MEs_col)
    annotation_row = data.frame(Module=module_name, row.names = module_name)
    module_name_without_me = substring(module_name,3)
    names(module_name_without_me) = module_name
    #print(as.data.frame(modTraitCor))
    sp_pheatmap(data=as.data.frame(modTraitCor), annotation_row = annotation_row,
                display_numbers = textMatrix,
                manual_annotation_colors_sidebar = list(Module=module_name_without_me),
                cluster_rows = T, cluster_cols = T,
                xtics_angle=45)

    if (!sp.is.null(saveplot)) {
      dev.off()
    }


    modTraitCorMelt = as.data.frame(modTraitCor)

    write.table(
      data.frame(ID = rownames(modTraitCor), modTraitCorMelt),
      file = paste0(prefix, ".module_trait_correlation.xls"),
      sep = "\t",
      quote = F, row.names=F
    )
    modTraitCorMelt$ID = rownames(modTraitCor)
    modTraitCorMelt = reshape2::melt(modTraitCorMelt)
    colnames(modTraitCorMelt) <-
      c("Module", "Trait", "PersonCorrelationValue")
    modTraitPMelt = as.data.frame(modTraitP)
    write.table(
      data.frame(ID = rownames(modTraitPMelt), modTraitPMelt),
      file = paste0(prefix, ".module_trait_correlationPvalue.xls"),
      sep = "\t",
      quote = F,
      row.names = F
    )
    modTraitPMelt$ID = rownames(modTraitP)
    modTraitPMelt = reshape2::melt(modTraitPMelt)
    colnames(modTraitPMelt) <- c("Module", "Trait", "Pvalue")
    #modTraitCorP = cbind(modTraitCorMelt, Pvalue=modTraitPMelt$Pvalue)
    modTraitCorP = merge(modTraitCorMelt, modTraitPMelt, by = c("Module", "Trait"))
    write.table(
      modTraitCorP,
      file = paste0(prefix, ".module_trait_correlationPvalueMelt.xls"),
      sep = "\t",
      quote = F,
      row.names = F
    )
    # checkAndInstallPackages(list(package1=c("ggmatrix","houyunhuang/ggmatrix")))
    #
    # modTraitCorP_ggmatrix <-
    #   ggplot(
    #     modTraitCorP,
    #     aes(
    #       x = Module,
    #       y = Trait,
    #       fill.upper = PersonCorrelationValue,
    #       fill.lower = Pvalue
    #     )
    #   ) +
    #   ggmatrix::geom_triangle() +
    #   ggmatrix::scale_fill_upper_gradientn(colours = up_color) +
    #   ggmatrix::scale_fill_lower_gradientn(colours = down_color) +
    #   theme(
    #     legend.position = "top",
    #     axis.text.x = element_text(
    #       angle = angle_x,
    #       hjust = 1,
    #       vjust = 0.3
    #     )
    #   ) + labs(x = "Module", y = "Trait")
    # ggsave(
    #   paste0(prefix, ".modTraitCorP_ggmatrix.pdf"),
    #   plot = modTraitCorP_ggmatrix,
    #   width = 15,
    #   height = 15,
    #   units = "cm"
    # )
    cat(sp_current_time(), "Finish plotting correlation heatmap among all traits and modules.\n")
  }


#' Plot gene module relationship and correlation with traits.
#'
#' @inheritParams WGCNA_cytoscape
#' @inheritParams WGCNA_saveModuleAndMe
#' @inheritParams WGCNA_coexprNetwork
#' @inheritParams WGCNA_MEs_traitCorrelationHeatmap
#'
#' @return A dataframe geneTraitCor
#' @export
#'
#' @examples
#'
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' traitData <- wgcnaL$traitData
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#' WGCNA_MEs_traitCorrelationHeatmap(MEs_col, traitData=traitData)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#' WGCNA_moduleTraitPlot(MEs_col, traitData=traitData)
#' geneTraitCor <- WGCNA_ModuleGeneTraitHeatmap(datExpr, traitData, net)
#'
#'
WGCNA_ModuleGeneTraitHeatmap <-
  function(datExpr,
           traitData=NULL,
           net,
           corType = "bicor",
           prefix = "ehbio",
           saveplot = NULL,
           ...) {
    # 计算性状与基因的相关性矩阵

    ## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵
    cat(sp_current_time(), "Plot module-gene-trait correlation heatmap.\n")
    if (sp.is.null(traitData)) {
      cat(sp_current_time(), "No module-gene-trait correlation heatmap plot since no trait data.\n")
      return(1)
    }else{
      traitData <- WGCNA_filterTrait(datExpr, traitData)
    }

    robustY = ifelse(corType == "pearson", T, F)
    if (corType == "pearson") {
      geneTraitCor = as.data.frame(cor(datExpr, traitData, use = "p"))
      nSamples <- nrow(traitData)
      geneTraitP = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneTraitCor), nSamples))
    } else {
      geneTraitCorA = WGCNA::bicorAndPvalue(datExpr, traitData, robustY = robustY)
      geneTraitCor = as.data.frame(geneTraitCorA$bicor)
      geneTraitP   = as.data.frame(geneTraitCorA$p)
    }

    geneTraitCorMelt = as.data.frame(geneTraitCor)
    write.table(
      geneTraitCorMelt,
      file = paste0(prefix, ".gene_trait_correlation.xls"),
      sep = "\t",
      quote = F
    )
    geneTraitCorMelt$ID = rownames(geneTraitCor)
    geneTraitCorMelt = reshape2::melt(geneTraitCorMelt)
    colnames(geneTraitCorMelt) <-
      c("Gene", "Trait", "PersonCorrelationValue")
    geneTraitPMelt = as.data.frame(geneTraitP)
    write.table(
      geneTraitPMelt,
      file = paste0(prefix, ".gene_trait_correlationPvalue.xls"),
      sep = "\t",
      quote = F
    )

    geneTraitPMelt$ID = rownames(geneTraitP)
    geneTraitPMelt = reshape2::melt(geneTraitPMelt)
    colnames(geneTraitPMelt) <- c("Gene", "Trait", "Pvalue")
    #geneTraitCorP = cbind(geneTraitCorMelt, Pvalue=geneTraitPMelt$Pvalue)
    geneTraitCorP = merge(geneTraitCorMelt, geneTraitPMelt, by = c("Gene", "Trait"))
    write.table(
      geneTraitCorP,
      file = paste0(prefix, ".gene_trait_correlationPvalueMelt.xls"),
      sep = "\t",
      quote = F,
      row.names = F
    )

    #geneTraitCorP = merge(geneTraitCorMelt, geneTraitPMelt, by=c("Gene","Trait"))
    write.table(
      geneTraitCorP[geneTraitCorP$Pvalue < 0.01, ],
      file = paste0(prefix, ".gene_trait_correlationPvalueMelt.p0.01.xls"),
      sep = "\t",
      quote = F,
      row.names = F
    )

    #plot_me_trat <- cbind(dynamicColors,moduleColors,geneTraitCor)
	#print(geneTraitCor)
    geneTraitCorColor <- WGCNA::numbers2colors(geneTraitCor)

    if (!sp.is.null(saveplot)) {
      base_plot_save(saveplot, bg = "white", ...)
    }

    moduleLabels = net$colors
    moduleColors = WGCNA::labels2colors(moduleLabels)
    dynamicColors <- WGCNA::labels2colors(net$unmergedColors)

    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]],
      cbind(dynamicColors, moduleColors, geneTraitCorColor),
      c("Dynamic Tree Cut", "Module colors", colnames(geneTraitCor)),
      dendroLabels = FALSE,
      hang = 0.5,
      addGuide = TRUE,
      guideHang = 0.05
    )

    if (!sp.is.null(saveplot)) {
      dev.off()
    }
    cat(sp_current_time(), "Finish plotting module-gene-trait correlation heatmap.\n")
    invisible(geneTraitCor)
  }


#' Genes correlated with both traits an modules.
#'
#' @inheritParams WGCNA_moduleTraitPlot
#' @param geneTraitCor A dataframe generated by \code{\link{WGCNA_ModuleGeneTraitHeatmap}}
#' @inheritParams WGCNA_saveModuleAndMe
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF(nSample=30, nGrp=3, sd=5)
#' datExpr <- WGCNA_dataFilter(df)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' WGCNA_saveModuleAndMe(net, datExpr)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#'
#' #2
#' exprMat <- "test.file"
#' wgcnaL <- WGCNA_readindata(exprMat)
#'
#' traitData <- 'trait.file'
#' wgcnaL <- WGCNA_readindata(exprMat, traitData)
#' datExpr <- wgcnaL$datExpr
#' traitData <- wgcnaL$traitData
#' WGCNA_dataCheck(datExpr)
#' datExpr <- WGCNA_dataFilter(datExpr)
#' datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr)
#' # datExpr <- WGCNA_sampleClusterDetectOutlier(datExpr, traitColors=wgcnaL$traitColors)
#' power <- WGCNA_softpower(datExpr)
#' net <- WGCNA_coexprNetwork(datExpr, power)
#' MEs_col <- WGCNA_saveModuleAndMe(net, datExpr)
#' WGCNA_MEs_traitCorrelationHeatmap(MEs_col, traitData=traitData)
#' cyt <- WGCNA_cytoscape(net, power, datExpr)
#' hubgene <- WGCNA_hubgene(cyt)
#' WGCNA_moduleTraitPlot(MEs_col, traitData=traitData)
#' geneTraitCor <- WGCNA_ModuleGeneTraitHeatmap(datExpr, traitData, net)
#' WGCNA_GeneModuleTraitCoorelation(datExpr, MEs_col, geneTraitCor, traitData, net)
#'
WGCNA_GeneModuleTraitCoorelation <-
  function(datExpr,
           MEs_col,
           geneTraitCor,
           traitData,
           net,
           corType = "bicor",
           prefix = "ehbio",
           ...) {
    ### 计算模块与基因的相关性矩阵

    if (sp.is.null(traitData)) {
      cat(sp_current_time(), "No module-gene-trait correlation plot since no trait data.\n")
      return(1)
    }else{
      traitData <- WGCNA_filterTrait(datExpr, traitData)
    }
    cat(sp_current_time(), "Plot interest genes for each module-trait combination group.\n")
    robustY = ifelse(corType == "pearson", T, F)
    if (corType == "pearson") {
      geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
      nSamples <- nrow(traitData)
      MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
    } else {
      # 关联样品性状的二元变量时，设置
      geneModuleMembershipA = WGCNA::bicorAndPvalue(datExpr, MEs_col, robustY =
                                                      robustY)
      geneModuleMembership = geneModuleMembershipA$bicor
      MMPvalue   = geneModuleMembershipA$p
    }


    geneModuleMembershipMelt = as.data.frame(geneModuleMembership)
    write.table(
      geneModuleMembershipMelt,
      file = paste0(prefix, ".gene_module_correlation.xls"),
      sep = "\t",
      quote = F
    )
    geneModuleMembershipMelt$ID = rownames(geneModuleMembership)
    geneModuleMembershipMelt = reshape2::melt(geneModuleMembershipMelt)
    colnames(geneModuleMembershipMelt) <-
      c("Gene", "Module", "PersonCorrelationValue")
    # geneTraitPMelt = as.data.frame(geneTraitP)
    # write.table(geneTraitPMelt,file=paste0(prefix,".gene_module_correlationMelt.xls"),
    #             sep="\t",quote=F)
    #
    # MMPvalueMelt$ID = rownames(MMPvalue)
    # MMPvalueMelt = reshape2::melt(MMPvalueMelt)
    # colnames(MMPvalueMelt) <- c("Gene","Module","Pvalue")
    # #geneModuleMembershipP = cbind(geneModuleMembershipMelt, Pvalue=MMPvalueMelt$Pvalue)
    # geneModuleMembershipP = merge(geneModuleMembershipMelt, MMPvalueMelt, by=c("Gene","Trait"))
    # write.table(geneModuleMembershipP,
    #             file=paste0(prefix,".gene_module_correlationPvalueMelt.xls"),
    #             sep="\t",quote=F,row.names=F)


    modNames = substring(colnames(MEs_col), 3)
    # 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
    #module = "red"
    #pheno = "Insulin_ug_l"

    phenoName = colnames(traitData)
    moduleLabels = net$colors
    moduleColors = WGCNA::labels2colors(moduleLabels)


    for (module in modNames) {
      if (module == "grey") {
        next
      }
      for (pheno in phenoName) {
        # 获取关注的列
        module_column = match(module, modNames)
        pheno_column = match(pheno, colnames(traitData))
        # 获取模块内的基因
        moduleGenes = moduleColors == module

        file <-
          paste(prefix,
                'gene_module_trait_cor',
                module,
                pheno,
                'xls',
                sep = ".")
        gene_trait_module_cor <-
          data.frame(ID=rownames(geneModuleMembership[moduleGenes,  module_column,drop=F]),
            geneModuleMembership = geneModuleMembership[moduleGenes,  module_column],
                geneTraitCor = geneTraitCor[moduleGenes,  pheno_column])
        #gene_trait_module_cor = data.frame(ID = rownames(gene_trait_module_cor), gene_trait_module_cor)
        write.table(
          gene_trait_module_cor,
          file = file,
          sep = "\t",
          quote = F,
          row.names = F
        )


        base_plot_save(paste(
          prefix,
          'gene_module_trait_cor',
          module,
          pheno,
          'pdf',
          sep = "."
        ),
        bg = "white",
        ...)

        sizeGrWindow(7, 7)
        par(mfrow = c(1, 1))
        # 与性状高度相关的基因，也是与性状相关的模型的关键基因
        WGCNA::verboseScatterplot(
          abs(geneModuleMembership[moduleGenes, module_column]),
          geneTraitCor[moduleGenes, pheno_column],
          #abs(geneTraitCor[moduleGenes, pheno_column]),
          xlab = paste("Module Membership in", module, "module"),
          ylab = paste("Gene significance for", pheno),
          main = paste("Module membership vs. gene significance\n"),
          cex.main = 1.2,
          cex.lab = 1.2,
          cex.axis = 1.2,
          col = module
        )

        dev.off()

        cat(sp_current_time(), paste0("Finish plotting interest genes for module (", module,
                                      ") - trait (", pheno, ") combination group.\n"))
      }
    }
  }



#' WGCNA onestep
#'
#' @inheritParams WGCNA_readindata
#' @inheritParams dataFilter
#' @inheritParams WGCNA_sampleClusterDetectOutlier
#' @inheritParams WGCNA_softpower
#' @inheritParams WGCNA_coexprNetwork
#' @inheritParams WGCNA_cytoscape
#' @inheritParams WGCNA_moduleTraitPlot
#' @param os_system Default the program will detect system type to choose which multiple thread function will be used.
#' `enableWGCNAThreads` is recommended, but only work in some linux os.
#' \code{\link[WGCNA]{allowWGCNAThreads}} is not recommended if `enableWGCNAThreads` works.
#' However for max and windows os, this is the only one can be used.
#' Even for some linux system, using `enableWGCNAThreads` will make programs
#' stuck in \code{\link[WGCNA]{pickSoftThreshold}} step.
#' So if stucked, supply any string other than `linux` to enable the usages
#' of \code{\link[WGCNA]{allowWGCNAThreads}}.
#' @param power_min For some data type, default selected power is a small number. Mostly this is due to unnormalized expression value, batch effects or small amount of total samples. When this happens, we may want to assign a power as 6 or other common numbers for downstream analysis. Here is where to specify it. Be careful to use this parameter unless you know what you are doing.
#'
#' @return net
#' @export
#'
#' @examples
#'
#' exprMat <- "test.file"
#'
#' traitData <- 'trait.file'
#'
#' WGCNA_onestep(exprMat, traitData)
#'
WGCNA_onestep <-
  function(exprMat,
           traitData = NULL,
           categoricalTrait = NULL,
           prefix = "ehbio",
           corType = "bicor",
           networkType = "signed",
           maxPower = NULL,
           maxBlockSize = NULL,
           top_mad_n = 0.75,
           rmVarZero = T,
           minimal_mad = NULL,
           thresholdZ.k = -2.5,
           TOM_plot = NULL,
           top_hub_n = 20,
           removeOutlier = F,
           RsquaredCut = 0.85,
           minModuleSize = NULL,
           mergeCutHeight = 0.2,
           numericLabels = TRUE,
           pamRespectsDendro = FALSE,
           saveTOMs = TRUE,
           maxPOutliers = NULL,
           loadTOM = TRUE,
           TOMDenom = "min",
           deepSplit = 1,
           stabilityCriterion = "Individual fraction",
           verbose = 0,
           os_system = NULL,
           randomSeed = 11521,
           dynamicCutPlot = TRUE,
           power_min = NULL,
           up_color = c("red", "white", "blue"),
           down_color = c("green", "white"),
           ...) {
    options(stringsAsFactors = FALSE)
    options(warn = -1)
    if (sp.is.null(os_system)) {
      os_system = Sys.info()['sysname']
    }

    #if (os_system == "Linux") {
    #  # 打开多线程
    #  enableWGCNAThreads()
    #} else {
    #  # if mac
    #  allowWGCNAThreads()
    #}


    wgcnaL <- WGCNA_readindata(exprMat, traitData = traitData,
                               categoricalTrait = categoricalTrait)

    datExpr <- wgcnaL$datExpr


    WGCNA_dataCheck(datExpr,
                    saveplot = paste0(prefix, ".WGCNA_dataCheck.pdf"),
                    width = 20)


    datExpr <-
      WGCNA_dataFilter(
        datExpr,
        minimal_mad = minimal_mad,
        top_mad_n = top_mad_n,
        rmVarZero = rmVarZero
      )


    datExpr <-
      WGCNA_sampleClusterDetectOutlier(
        datExpr,
        traitColors = wgcnaL$traitColors,
        thresholdZ.k = thresholdZ.k,
        removeOutlier = removeOutlier,
        saveplot = paste0(prefix, ".WGCNA_sampleClusterDetectOutlier.pdf")
      )



    power <-
      WGCNA_softpower(
        datExpr,
        saveplot = paste0(prefix, ".WGCNA_softpower.pdf"),
        networkType = networkType,
        maxPower = maxPower,
        RsquaredCut = RsquaredCut
      )

    power <- power$power

    if (!sp.is.null(power_min) && (power < power_min)) {
      power = power_min
    }

    net <- WGCNA_coexprNetwork(
      datExpr,
      power,
      saveplot = paste0(prefix, ".WGCNA_module_generation_plot.pdf"),
      maxBlockSize = maxBlockSize,
      minModuleSize = minModuleSize,
      networkType = networkType,
      mergeCutHeight = mergeCutHeight,
      numericLabels = numericLabels,
      pamRespectsDendro = pamRespectsDendro,
      saveTOMs = saveTOMs,
      corType = corType,
      maxPOutliers = maxPOutliers,
      loadTOM = loadTOM,
      TOMDenom = TOMDenom,
      deepSplit = deepSplit,
      stabilityCriterion = stabilityCriterion,
      saveTOMFileBase = paste0(prefix, ".blockwiseTOM"),
      verbose = verbose,
      randomSeed = randomSeed,
      dynamicCutPlot = dynamicCutPlot
    )

    MEs_col <- WGCNA_saveModuleAndMe(
      net,
      datExpr,
      prefix = prefix,
      saveplot = paste0(prefix, ".WGCNA_module_correlation_plot.pdf")
    )


    net$MEs_col <- MEs_col

    # WGCNA_MEs_traitCorrelationHeatmap(
    #   MEs_col,
    #   traitData = traitData,
    #   saveplot = paste0(prefix, ".WGCNA_moduletrait_correlation_plot.pdf")
    # )

    cyt <-
      WGCNA_cytoscape(net, power, datExpr, TOM_plot = TOM_plot, prefix = prefix)

    net$cyt <- cyt

    hubgene <- WGCNA_hubgene(cyt, top_hub_n = top_hub_n, prefix = prefix)

    net$hubgene <- hubgene

    if (!is.null(traitData)) {
      WGCNA_moduleTraitPlot(
        MEs_col,
        traitData = wgcnaL$traitData,
        saveplot = paste0(prefix, ".WGCNA_moduleTraitHeatmap.pdf"),
        corType = corType,
        prefix = prefix,
        ...
      )

      geneTraitCor <-
        WGCNA_ModuleGeneTraitHeatmap(
          datExpr,
          traitData = wgcnaL$traitData,
          net = net,
          prefix = prefix,
          saveplot = paste0(prefix, ".WGCNA_ModuleGeneTraitHeatmap.pdf")
        )

      net$geneTraitCor <- geneTraitCor

      WGCNA_GeneModuleTraitCoorelation(
        datExpr,
        MEs_col,
        geneTraitCor,
        traitData = wgcnaL$traitData,
        net,
        corType = corType,
        prefix = prefix
      )
    }
    invisible(net)
    cat(sp_current_time(), "Success.\n")
  }

### 计算邻接矩阵
# adjacency = adjacency(datExpr, power = power)
#
# ### 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
# TOM = TOMsimilarity(adjacency)
# dissTOM = 1-TOM
#
# ### 层级聚类计算基因之间的距离树
# geneTree = flashClust(as.dist(dissTOM), method = "average")
#
# ### 模块合并
# # We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30
# # Module identification using dynamic tree cut:
# dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#                             deepSplit = 4, pamRespectsDendro = FALSE,
#                             minClusterSize = minModuleSize)
# # Convert numeric lables into colors
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
#
# ### 通过计算模块的代表性模式和模块之间的定量相似性评估，合并表达图谱相似
# #的模块
# MEList = moduleEigengenes(datExpr, colors = dynamicColors)
# MEs = MEList$eigengenes
# # Calculate dissimilarity of module eigengenes
# MEDiss = 1-cor(MEs)
# # Cluster module eigengenes
# METree = flashClust(as.dist(MEDiss), method = "average")
# MEDissThres = 0.25
#
# # Call an automatic merging function
# merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# # The merged module colors
# mergedColors = merge$colors
# table(mergedColors)
# # Eigengenes of the new merged
#
# plotDendroAndColors(geneTree, cbind(dynamicColors,moduleColors),
#                     c("Dynamic Tree Cut", "Module colors"),
#                     dendroLabels = FALSE, hang = 0.5,
#                     addGuide = TRUE, guideHang = 0.05)

## 分步法完结
