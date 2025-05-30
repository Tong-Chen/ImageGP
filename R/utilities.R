


# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'
#   Generate DOC (Mac):              'Ctrl + Shift + Option + r'

#' Check and install given packages
#'
#' @param package A list containing names and install-names of each package.
#' (instll-names is only required for packages from github.)
#' Like list(package1=c("ggplot2")) or
#' list(packages1=c("ggplot2"), package2=c("ImageGP", "git_user/ImageGP"))
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' checkAndInstallPackages(list(package1=c("ggplot2")))
#'
#' checkAndInstallPackages(list(packages1=c("ggplot2"), package2=c("ImageGP", "git_user/ImageGP")))
#'
checkAndInstallPackages <-
  function(packageL, site = "https://mirrors.tuna.tsinghua.edu.cn/CRAN") {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager",
                       update = F,
                       site_repository = site)

    all_installed_packages = rownames(installed.packages())
    for (i in packageL) {
      package_name = i[1]
      package_install_name = i[1]
      if (length(i) == 2) {
        package_install_name = i[2]
      }
      if (!package_name %in% all_installed_packages)
        BiocManager::install(package_install_name,
                             update = F,
                             site_repository = site)
    }

    suppressPackageStartupMessages(library(package_name, character.only = TRUE))
  }


#' Extract same rows of two data.frame with same order.
#'
#' @param df1 Dataframe 1
#' @param df2 Dataframe 2
#' @param way `row-row (default)`: Extract same rows of two data.frame with same order.
#'            `col-row`: Extarct same columns in df1 with rows in df2 and with same order
#'
#' @return A list in format like list(df1=df1, df2=df2)
#' @export
#'
#' @examples
#'
#' NULL
#'
match_two_df <- function(df1, df2, way="row-row"){
  if (way == "row-row"){
    df1_rownameL <- rownames(df1)
    df2_rownameL <- rownames(df2)

    common_rownameL <- intersect(df1_rownameL, df2_rownameL)

    # 保证表达表样品与df2样品顺序和数目完全一致
    df1 <- df1[common_rownameL,,drop=F]
    df2 <- df2[common_rownameL,,drop=F]
  } else if (way == "col-row") {
    df1_colnameL <- colnames(df1)
    df2_rownameL <- rownames(df2)

    common_rownameL <- intersect(df1_colnameL, df2_rownameL)

    # 保证表达表样品与df2样品顺序和数目完全一致
    df1 <- df1[,common_rownameL,drop=F]
    df2 <- df2[common_rownameL,,drop=F]
  }
  return(list(df1=df1, df2=df2))
}

#' Get current time in strign format
#'
#' @param delim_left Default `[`.
#' @param delim_right Default `]`.
#'
#' @return A string
#' @export
#'
#' @examples
#'
#' sp_current_time()
#'
sp_current_time <- function(delim_left = '[',
                            delim_right = ']') {
  return(paste0(delim_left, Sys.time(), delim_right))
}

#' Determine the value to add befor log transform.
#'
#' @param data A numerical dataframe or a vector
#' @param ratio Minimum non-zero value would be used as add values. if `ratio` specified,
#' the detected minimum non-zero multiple ratio would be returned.
#'
#' @return A numericalvalue
#' @export
#'
#' @examples
#'
#' sp_determine_log_add(c(1,2,3))
#'
sp_determine_log_add <- function(data, ratio = 1) {
  min_value = min(min(data, na.rm=T), na.rm=T)
  if (min_value > 0) {
    return(0)
  } else if (min_value == 0) {
    min_value = min(min(data[data != 0],na.rm=T))
    return(min_value * ratio)
  } else{
    stop("Negative value is not allowed for log2 transform!")
  }
}

#' Check Null Object
#'
#' @param x `NULL` object or `'null'` string
#'
#' @return True when x is `NULL` or `"NULL"` (case insensitive for character type)
#' @export
#'
#' @examples
#'
#' sp.is.null('NULL')
#'
sp.is.null <- function(x) {
  if (length(x) > 1) {
    return(FALSE)
  }
  if (is.character(x)) {
    return(toupper(x) == 'NULL')
  } else{
    return(base::is.null(x))
  }
}


#' Transfer color string to vector
#'
#' @param x A string
#' @param pattern delimiter of sub-strings
#'
#' @return A vector
#' @export
#'
#' @examples
#'
#' sp_string2vector('red, blue,white')
#'
sp_string2vector <- function(x, pattern = ",") {

  if(sp.is.null(x)){
    return(x)
  }
  if (requireNamespace("stringr", quietly = TRUE)) {
    str2v <- stringr::str_trim(stringr::str_split(x, pattern, simplify = T))
  } else {
    str2v <- trimws(unlist(strsplit(x, split = pattern)))
  }
  if(numCheck(str2v)){
  	str2v <- mixedToFloat(str2v)
  }
  return(str2v)
}

#' Read in data
#'
#' @inheritParams  utils::read.table
#' @param renameDuplicateRowNames If TRUE, the function will transfer first column
#' as row names (with duplicates numbered)
#' @param ... Other parameters given to \code{read.table}
#'
#' @return data.frame
#' @export
#'
#' @examples
#'
#' # Not run
#' sp_readTable("a.txt")
#'
sp_readTable <-
  function(file,
           sep = "\t",
           row.names = NULL,
           header = T,
           quote = "",
           comment = "",
           check.names = F,
           renameDuplicateRowNames = F,
           stringsAsFactors = T,
           ...) {
    if (renameDuplicateRowNames) {
      data <- read.table(
        file,
        sep = sep,
        row.names = NULL,
        header = header,
        quote = quote,
        comment = comment,
        check.names = check.names,
        stringsAsFactors = stringsAsFactors,
        ...
      )
      if (!is.numeric(as.vector(data[, 1]))){
        rownames_data <- make.unique(as.vector(data[, 1]))
      } else {
        rownames_data <- as.vector(data[, 1])
      }

      data <- data[,-1, drop = F]
      rownames(data) <- rownames_data
    } else {
      data <-
        read.table(
          file,
          sep = sep,
          row.names = row.names,
          header = header,
          quote = quote,
          comment = comment,
          check.names = check.names,
          stringsAsFactors = stringsAsFactors,
          ...
        )
    }
    invisible(data)
  }

#' Write dataframe to file with names of first column filled.
#'
#' @param df A dataframe
#' @param file Filename
#' @param keep_rownames Default TRUE meaning output rownames as the first column
#' with column name is \code{ID}. If FALSE, ignore rownames.
#' @inheritParams utils::write.table
#'
#' @return NA
#' @export
#'
#' @examples
#'
#' # Not run
#' sp_writeTable(df, "a.txt")
#'
sp_writeTable <- function(df,
                          file = '',
                          keep_rownames = T,
                          col.names=T) {
  if (keep_rownames) {
    write.table(
      data.frame(ID = rownames(df), df),
      file = file,
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = col.names
    )
  } else {
    write.table(
      df,
      file = file,
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = col.names
    )
  }
}

#' Generate gene expression table or otu abundance table with given samle information for test.
#'
#' @param type Generate gene expression or OTU abundance. Only affect rownames.
#' @param mean Mean value of abundance given to \code{\link{rnorm}}.
#' @param sd Standard deviations given to \code{\link{rnorm}}.
#' @param nGene Number of genes or OTUs.
#' @param nGrp Number of sample groups.
#' @param nSample Number of sample replications for each group.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#'
generateAbundanceDF <-
  function(type = "Gene",
           mean = 20,
           nGene = 15,
           nGrp = 2,
           nSample = 3) {
    df <-
      as.data.frame(matrix(rnorm(nGene * nGrp * nSample, mean = mean), nrow =
                             nGene))
    colnames(df) <-
      paste("Samp", paste(rep(LETTERS[1:nGrp], each = nSample), rep(1:nSample, nGrp), sep =
                            "_"), sep = "_")
    rownames(df) <- paste(type, letters[1:nGene], sep = "_")
    return(df)
  }



#' Get ordered column correlation matrix from input dataframe. Normally used
#' to do sample corealtion of gene expression or OTU abundance matrix.
#'
#' @param mat A dataframe.
#' @param method Type of correlation coefficient given to \code{\link{cor}}.
#' Default "pearson".
#' @param digits Number of decimial digits (given to \code{\link{round}}) to keep (default 4).
#' @param cor_file Save ordered correlation matrix to given file name.
#'
#' @return A list containing ordered column correlation matrix and hcluster result.
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' Matrix2colCorrelation(df)
#'
Matrix2colCorrelation <-
  function(mat,
           method = "pearson",
           digits = 4,
           cor_file = NULL) {
    pearson_cor <-
      round(as.matrix(cor(mat, method = method)), digits = digits)
    hc <- amap::hcluster(t(mat), method = method)
    pearson_cor <- pearson_cor[hc$order, hc$order]
    if (!is.null(file)) {
      pearson_cor_output = data.frame(id = rownames(pearson_cor), pearson_cor)
      write.table(
        pearson_cor_output,
        file = cor_file,
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = T
      )
    }
    return(list(pearson_cor = pearson_cor, hc = hc))
  }



#' Get lower triangle of the correlation matrix (from web)
#'
#' @param cormat A data frame
#'
#' @return A data frame
#' @export
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' df_cor = Matrix2colCorrelation(df)
#' get_lower_tri(df_cor)
#'
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}


#' Get upper triangle of the correlation matrix (from web)
#'
#' @param cormat A data frame
#'
#' @return A data fram
#' @export
#'
#'
#' @examples
#'
#' df = generateAbundanceDF()
#' df_cor = Matrix2colCorrelation(df)
#' get_upper_tri(df_cor)
#'
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}



# options(scipen=999)

#' Check if given string or vector is all numeric
#'
#' @param x A string or a vector
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#'
#' numCheck(3)
#'
#' numCheck("-1/3")
#'
#' numCheck(c("1","0.2","1/3","-1"))
#'
numCheck <- function(x) {
  # Function get from https://stackoverflow.com/questions/10674992/convert-a-character-vector-of-mixed-numbers-fractions-and-integers-to-numeric?rq=1
  # With little modifications
  is.numeric2 <- is.numeric(x)
  if(is.numeric2){
    return(is.numeric2)
  }
  is_na <- is.na(x)
  x <- sapply(x, as.character)
  x <- trimws(x)
  is.integer  <- grepl("^-?\\d+$", x)
  is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
  is.float <- grepl("^-?\\d+\\.\\d+$", x)
  is.percent <- grepl("[0-9.]+%$", x)
  is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
  return(all(
      is.integer | is.fraction | is.float | is.mixed | is.percent | is_na
  ))
}

#' Transfer numeric string to numeric.
#'
#' @param x A string or a vector
#'
#' @return A number or a numeric vector
#' @export
#'
#' @examples
#'
#' mixedToFloat(3)
#'
#' mixedToFloat("-1/3")
#'
#' mixedToFloat(c("1","0.2","1/3","-1"))
#'
mixedToFloat <- function(x) {
  is_na <- is.na(x)
  x <- sapply(x, as.character)
  x <- trimws(x)
  is.integer  <- grepl("^-?\\d+$", x)
  is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
  is.float <- grepl("^-?\\d+\\.\\d+$", x)
  is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
  is.percent <- grepl("[0-9.]+%$", x)
  stopifnot(all(is.integer |
                  is.fraction | is.float | is.mixed | is.percent | is_na))

  numbers <- strsplit(x, "[ /%]")

  ifelse(is_na,
         NA,
         ifelse(
           is.integer,
           as.numeric(sapply(numbers, `[`, 1)),
           ifelse(
             is.percent,
             as.numeric(sapply(numbers, `[`, 1)) / 100,
             ifelse(
               is.float,
               as.numeric(sapply(numbers, `[`, 1)),
               ifelse(
                 is.fraction,
                 as.numeric(sapply(numbers, `[`, 1)) /
                   as.numeric(sapply(numbers, `[`, 2)),
                 as.numeric(sapply(numbers, `[`, 1)) +
                   as.numeric(sapply(numbers, `[`, 2)) /
                   as.numeric(sapply(numbers, `[`, 3))
               )
             )
           )
         ))

}

#mixedToFloat(c('1 1/2', '2 3/4', '2/3', '11 1/4', '1'))

#' GGSCI object for json
#'
#' @return A color vector
#' @export
#'
#' @examples
#'
#' ggsci_to_json()
#'
ggsci_to_json <- function(){
  library(jsonlite)

  json_list <- list()

  json_data <- lapply(names(ggsci_db), function(group) {
    #cat(group)
    lapply(names(ggsci_db[[group]]), function(subgroup){
      #cat(group)
      #cat(subgroup)
      colors <- unname(ggsci_db[[group]][[subgroup]])
      short_name <- paste("ggsci", group, subgroup, sep = "_")
      long_name <- paste("ggsci", group, subgroup, sep = "_")
      json_list[[length(json_list) + 1]] <<- list(color = colors,
                                                  short_name = unbox(short_name),
                                                  long_name = unbox(long_name),
                                                  type = unbox("string"))
    })

  })

  # Convert the list to JSON
  json_output <- toJSON(json_list, pretty = TRUE)

  # Print the JSON
  # cat(json_output)

  writeLines(json_output,"ggsci_colors.json")
}




#' Generate color code
#'
#' @param color Colors like c('red', 'blue', '#6181BD') or
#' a RColorBrewer color set like  "BrBG"     "PiYG"     "PRGn"     "PuOr"
#' "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"  "Spectral" "Accent"
#' "Dark2"    "Paired"   "Pastel1"  "Pastel2"  "Set1"
#' "Set2"    "Set3"     "Blues"    "BuGn"     "BuPu"
#' "GnBu"     "Greens"   "Greys"    "Oranges" "OrRd"     "PuBu"
#' "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"
#' "YlGn"    "YlGnBu"   "YlOrBr"   "YlOrRd"
#' (check http://www.sthda.com/english/wiki/colors-in-r for more).
#'
#' @param number Number of colors to return.
#'
#' @param alpha Generate an alpha transparency values for return colors. 0 means fully transparent and 1 means opaque. Default 1.
#'
#' @param reverseColorList Get the reverse of generated color list.
#' @return A color vector
#' @export
#'
#' @examples
#'
#' generate_color_list('red', 5)
#' generate_color_list(c('green', 'red'), 5)
#' generate_color_list(c('green', 'red'), 5, alpha=0.5)
#' generate_color_list("Set3", 5)
#'



generate_color_list <- function(color, number, alpha = 1, constantColor=F, reverseColorList=F) {
  color = color[color!="None" & color !=""]
  color_len = length(color)

  if (color_len == 1) {
    if sp.is.null(color){
      color = "Set2"
    }
    if (grepl("^ggsci", color)) {
      # color = "ggsci_npg_nrc"
      colorV = sp_string2vector(color, "_")
      colorSet = unname(ggsci_db[[colorV[2]]][[colorV[3]]])
      maxcolor = length(colorSet)
      if (number <= maxcolor) {
        colorL <- colorSet[1:number]
      } else {
        colorL <-
          colorRampPalette(colorSet)(number)
      }
    } else {
      brewer = rownames(RColorBrewer::brewer.pal.info)
      if (color %in% brewer) {
        maxcolor = RColorBrewer::brewer.pal.info[color, ]$maxcolors
        if (number <= maxcolor) {
          mincolor = 3
          if (number < mincolor) {
            colorL <- RColorBrewer::brewer.pal(mincolor, color)[1:number]
          } else {
            colorL <- RColorBrewer::brewer.pal(number, color)
          }
        } else {
          colorL <-
            colorRampPalette(RColorBrewer::brewer.pal(maxcolor, color))(number)
        }
      } else{
        if(constantColor){
          colorL <- c(color, rep("gray23",number-1))
        } else {
          colorL <- rep(color, number)
        }
      }
    }

  } else if (color_len == number) {
    colorL = color
    # return(colorL)
  } else{
    if(constantColor) {
      if ((number - color_len) < 0) {
        colorL <- color[1:number]
      } else{
        colorL <- c(color, rep("gray23", number - color_len))
      }
    } else {
      colorL = colorRampPalette(color, alpha=T)(number)
    }
    # return(colorL)
  }

  if (reverseColorList){
    colorL = rev(colorL)
  }

  if (alpha==1){
    return(colorL)
  }

  return(rgb(
    t(col2rgb(colorL)),
    alpha = alpha * 255,
    maxColorValue = 255
  ))
}

#' Transfer one column of data.
#'
#' @param data A data matrix
#' @param variable One column name of data matrix
#' @param y_add A number to add if log scale is used.
#' Default 0 meaning the minimum non-zero value would be used.
#' @param yaxis_scale_mode Give the following `scale_y_log10()`,
#' `coord_trans(y="log10")`, or other legal command for ggplot2 or
#' simply `log2` to set the scale way.
#'
#' @return A data frame
#' @export
#'
#' @examples
#'
#' data <- data.frame(A=letters[1:4], B=letters[1:4])
#' data
#' data = sp_transfer_variable(data,'A', "log2")


sp_transfer_one_column <- function(data, variable, yaxis_scale_mode=NULL, y_add=0){
  if(numCheck(data[[variable]])){
    if (!is.numeric(data[[variable]])) {
      data[[variable]] <- mixedToFloat(data[[variable]])
    }
  } else {
    stop(paste(variable,"column is not numerical column."))
  }
  # print(y_add)
  # Give the minimum non-zero value to add to avoid log2(0)
  if (y_add == 0) {
    y_add = sp_determine_log_add(data[[variable]])
    # print(paste("153", y_add))
  }
  # print("155")
  # print(data[[yvariable]])
  data[[variable]] <- data[[variable]] + y_add
  if (yaxis_scale_mode == "log2") {
    data[[variable]] <- log2(data[[variable]])
  } else if (yaxis_scale_mode == "log10") {
    data[[variable]] <- log10(data[[variable]])
  }
  return(data)
}



#' Set factor order of given variable. If `variable_order` is supplied, only
#' factors in `variable_order` will be kept and re-factored. Other variables
#' would be depleted.
#'
#' @param data A data matrix
#' @param variable One column name of data matrix
#' @param variable_order Expected order of `data[[variable]]`.
#' @param order_data_frame_by_this_variable_order Return ordered dataframe by this order.
#' Please remember that only keep the last order if applying multiple order operation.
#' @param filter_unexist_factor Filter un-exist factors.
#' @param rename_levels Rename old levels to new levels. Default False.
#'
#' @return A data frame
#' @export
#'
#' @examples
#'
#' data <- data.frame(A=letters[1:4], B=letters[1:4])
#' data
#' data = sp_set_factor_order(data,'A')
#' data$A
#' data = sp_set_factor_order(data,'B',c('c','d','b','a'))
#' data$B
#' data = sp_set_factor_order(data,'B',c('c','d','a'))
#' data$B
#'
sp_set_factor_order <-
  function(data, variable, variable_order = NULL, order_data_frame_by_this_variable_order=F,
           filter_unexist_factor=T, rename_levels=F) {
    if (!variable %in% colnames(data)){
      stop(paste(variable,'must be one of column names of data!'))
    }
    if(numCheck(data[[variable]])){
      if (!is.numeric(data[[variable]])) {
        data[[variable]] <- mixedToFloat(data[[variable]])
      }
      # add to filter by numeric values
      if (!sp.is.null(variable_order) && length(variable_order)==2){
        variable_order <- sort(mixedToFloat(variable_order))
        data <- data[data[[variable]]>=variable_order[1] & data[[variable]]<=variable_order[2], ,drop=F]
        if(nrow(data)==0){
          stop(paste0("No data avaiable after filtering by column <", variable, "> with <",
                      paste(variable_order, collapse=","),">"))
        }
      }
    } else {
      if (!sp.is.null(variable_order)) {
        if (filter_unexist_factor){
          data = data[data[[variable]] %in% variable_order, , drop = F]
          data = droplevels(data)
        }
        if (rename_levels){
          levels(data[[variable]]) <- variable_order
        } else {
          data[[variable]] <-
            factor(data[[variable]], levels = variable_order, ordered = T)
        }

      } else {
        data[[variable]] <- factor(data[[variable]],
                                   levels = unique(data[[variable]]), ordered = T)
      }
    }
    if (order_data_frame_by_this_variable_order){
      data = data[rev(order(data[[variable]])),,drop=F]
    }
    invisible(data)
  }

#' Add manual color assignment for both categorical and numerical variable
#'
#' @param p A ggplot2 object
#' @param data Data matrix used for the ggplot2 object `p`
#' @param color_variable Name of columns for color assignment
#' @param manual_color_vector Manually set colors for each geom.
#' Default NULL, meaning using ggplot2 default.
#' Colors like c('red', 'blue', '#6181BD') (number of colors not matter) or
#' a RColorBrewer color set like  "BrBG"     "PiYG"     "PRGn"     "PuOr"
#' "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"  "Spectral" "Accent"
#' "Dark2"    "Paired"   "Pastel1"  "Pastel2"  "Set1"
#' "Set2"    "Set3"     "Blues"    "BuGn"     "BuPu"
#' "GnBu"     "Greens"   "Greys"    "Oranges" "OrRd"     "PuBu"
#' "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"
#' "YlGn"    "YlGnBu"   "YlOrBr"   "YlOrRd"
#' (check http://www.sthda.com/english/wiki/colors-in-r for more).
#' @param alpha Color transparency (0-1). 0: opaque; 1: transparent.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' p <- sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)
#'
#' ## End(Not run)
#'
sp_manual_color_ggplot2 <-
  function (p,
            data,
            color_variable,
            manual_color_vector = NULL,
            alpha = 1) {
    if (!sp.is.null(manual_color_vector)) {
      if (is.numeric(data[[color_variable]])) {
        color_v <- generate_color_list(manual_color_vector, 10, alpha = alpha)
        p <-
          p + scale_color_gradientn(colors = color_v)
      } else {
        color_v <-
          generate_color_list(manual_color_vector, length(unique(data[[color_variable]])),
                              alpha = alpha)
        p <- p + scale_color_manual(values = color_v)
      }
    }
    p
  }

#' Add manual fill-color assignment for both categorical and numerical variable
#'
#' @param p A ggplot2 object
#' @param data Data matrix used for the ggplot2 object `p`
#' @param color_variable Name of columns for color assignment
#' @inheritParams sp_manual_color_ggplot2
#' @param alpha Transparency
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' p <- sp_manual_fill_ggplot2(p, data, color_variable, manual_color_vector)
#'
#' ## End(Not run)
#'
sp_manual_fill_ggplot2 <-
  function (p,
            data,
            color_variable,
            manual_color_vector = NULL,
            alpha = 1) {
    if (!sp.is.null(manual_color_vector)) {
      if (is.numeric(data[[color_variable]])) {
        color_v <- generate_color_list(manual_color_vector, 10, alpha = alpha)
        print(color_v)
        p <-
          p + scale_fill_gradientn(colors = color_v)
      } else {

        color_v <-
          generate_color_list(manual_color_vector, length(unique(data[[color_variable]])),
                              alpha = alpha)
        print(color_v)
        p <- p + scale_fill_manual(values = color_v)
      }
    }
    p
  }

#' Add hline or vline for ggplot2 object
#'
#' @param p A ggplot2 object
#' @param custom_vline_x_position A vector of coordinates for vertical lines.
#' @param custom_vline_anno Annotation text for each vertical line.
#' @param custom_hline_y_position A vector of coordinates for horizontal lines.
#' @param custom_hline_anno Annotation text for each horizontal line.
#' @inheritParams ggplot2::geom_vline
#' @param ... Extra parameters given to `geom_vline` and `geom_hline`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' sp_ggplot_add_vline_hline(p)
#'
#' ## End(Not run)
#'
#'
sp_ggplot_add_vline_hline <- function(p,
                                      custom_vline_x_position = NULL,
                                      custom_vline_anno = NULL,
                                      custom_vline_anno_y_pos = NULL,
                                      custom_hline_y_position = NULL,
                                      custom_hline_anno = NULL,
                                      custom_hline_anno_x_pos = NULL,
                                      linetype = "dotted",
                                      size = 0.5,
                                      ...) {
  if (!sp.is.null(custom_vline_x_position) || !sp.is.null(custom_hline_y_position)){
    gb = ggplot_build(p)
  }

  if (!sp.is.null(custom_vline_x_position)) {
    p <- p + geom_vline(xintercept = custom_vline_x_position,
                        linetype = linetype,
                        linewidth = size,
                        ...)
    if (!is.null(custom_vline_anno)) {
      if (is.null(custom_vline_anno_y_pos)) {
        custom_vline_anno_y_pos = gb$layout$panel_params[[1]]$y.range[2]
      }
      p <-
        p + annotate(
          "text",
          x = custom_vline_x_position,
          y = custom_vline_anno_y_pos,
          label = custom_vline_anno,
          hjust = 0
        )
    }
  }


  if (!sp.is.null(custom_hline_y_position)) {
    p <- p + geom_hline(yintercept = custom_hline_y_position,
                        linetype = linetype,
                        linewidth = size,
                        ...)
    if (!is.null(custom_hline_anno)) {
      if (is.null(custom_hline_anno_x_pos)) {
        custom_hline_anno_x_pos = 0
      }
      #xmax = gb$layout$panel_params[[1]]$x.range[2]
      p <-
        p + annotate(
          "text",
          y = custom_hline_y_position,
          x = custom_hline_anno_x_pos,
          label = custom_hline_anno,
          vjust = 0,
          hjust = 0
        )
    }
  }
  return(p)
}

#' Facet ggplot2 object
#'
#' @param p A ggplot2 object
#' @param facet_variable Wrap plots by given column (one of column names should be specified).
#' This is used to put multiple plot in one picture.
#' @param facet_nrow 	The number of rows one want when `facet` is used. Default NULL.
#' @param facet_ncol The number of columns one want when `facet` is used. Default NULL.
#' @param facet_scales Paramter for scales for facet. Default `fixed` meaning each inner graph
#' use same scale (x,y range), `free` (variable x, y ranges for each sub-plot),
#' `free_x` (variable x ranges for each sub-plot), `free_y` (variable y ranges for each sub-plot).
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' sp_ggplot_facet(p, facet_variable)
#'
#' ## End(Not run)
#'

sp_ggplot_facet <- function(p, facet_variable=NULL, facet_ncol=NULL, facet_nrow=NULL, facet_scales="fixed"){
  p <- p + facet_wrap( ~  .data[[facet_variable]],
                       ncol = facet_ncol,
                       nrow = facet_nrow,
                       scales = facet_scales)
  return(p)
}


#' Used to read in long/wide format file or datafrmes. Wide format would be transferred to lonf fromat.
#'
#' @param data Data frame or data file (with header line, the first column will
#' not be treated as row names for long format matrix, tab seperated).
#' @param xvariable Name for x-axis variable.
#' @param melted `TRUE` for dealinig with long format matrix, the program will skip melt preprocess.
#' Default `FALSE` for dealing with wide format matrix.
#' @param ... Parameters given to \code{\link{dataFilter2}}.
#' @return a A long format dataframe
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#'
#' ## Not run:
#' sp_read_in_long_wide_matrix(data, xvariable, melted)
#'
#' ## End(Not run)
#'
sp_read_in_long_wide_matrix <- function(data, xvariable, melted, ...){
  wide_rownames = NULL
  wide_colnames = NULL
  if ("character" %in% class(data)) {
    if (!melted) {
      data <- sp_readTable(data, row.names = NULL)
      wide_rownames <- make.unique(as.vector(as.character(data[, 1])))
      data <- data[, -1, drop = F]
      rownames(data) <- wide_rownames
      wide_colnames <- colnames(data)

      if (all(apply(data, 2, numCheck), na.rm = T)) {
        rownames_data <- rownames(data)
        data <- as.data.frame(apply(data, 2, mixedToFloat))
        adata <- as.data.frame(data)
        rownames(data) <- rownames_data
      } else {
        stop(
          "For wide format data matrix, all elements except the first row and column must be numbers unless long format is used."
        )
      }
      # print(data)
      data <- dataFilter2(data, ...)
      # print(data)
      wide_rownames <- rownames(data)
      data[[xvariable]] <- wide_rownames
      data <- reshape2::melt(data, id.vars = xvariable)
    } else {
      data <- sp_readTable(data, row.names = NULL)
    }
  } else{
    if(!"data.frame" %in% class(data)) {
      stop("Unknown input format for `data` parameter.")
    }
    if (!melted) {
      wide_colnames <- colnames(data)
      data <- dataFilter2(data, ...)
      wide_rownames <- rownames(data)
      data[[xvariable]] <- wide_rownames
      data <- reshape2::melt(data, id.vars = xvariable)
    }
  }
  invisible(list(data=data, wide_rownames=wide_rownames, wide_colnames=wide_colnames))
}


#' Generate shape symbols for large number of groups for ggplot2
#'
#' @param data A data frame
#' @param shape_variable The variable treated as shape groups
#'
#' @return A vector contains all group symbols
#' @export
#'
#' @examples
#'
#' # Not run
#'
generate_shapes <- function(data, shape_variable){
  shape_level <- nlevels(data[[shape_variable]])
  if (shape_level < 15){
    shapes = (0:shape_level) %% 15
  } else{
    shapes = c(0:14,c((15:shape_level) %% 110 + 18))
  }
}

#' Use showtext to load fonts
#'
#' @param font_path Specify font type. Give a path for one font type file
#' like '/etc/fonts/Arial.ttf'
#' or 'HeiArial.ttc'(if in current directory), Default system default.
#'
#' @return font_name or null
#' @export
#'
#' @examples
#'
#' ## Not run:
#' sp_load_font(font_path="arial.tff")
#'
#' ## End(Not run)
#'
sp_load_font <- function(font_path){
  if (!sp.is.null(font_path)) {
    if (!requireNamespace("showtext", quietly = TRUE))
      install.packages("showtext", quite=T)
    library(showtext)
    showtext.auto(enable = TRUE)
    font_name = tools::file_path_sans_ext(basename(font_path))
    font.add(font_name, font_path)
    return(font_name)
  }
  return(NULL)
}



#' Change common layout of ggplot2 object
#'
#' @param p A ggplot2 object
#' @param xtics_angle Rotation angle for a-axis. Default 0.
#' @param legend.position Position of legend, accept top, bottom, left, right, none or c(0.8,0.8).
#' @param extra_ggplot2_cmd Extra ggplot2 commands (currently unsupported)
#' @param filename Output picture to given file.
#' @param title Title of picture.
#' @param x_label Xlab label.
#' @param y_label Ylab label.
#' @param coordinate_flip Flip cartesian coordinates so that horizontal becomes vertical, and vertical, horizontal. This is primarily useful for converting geoms and statistics which display y conditional on x, to x conditional on y.
#' @param width Picture width (units: cm)
#' @param height Picture height (units: cm)
#' @param zoom_split If both x and y is given, should each axis zoom be shown separately as well? Defaults to FALSE.
#' @param zoom_xlim Specific zoom ranges for x axis.
#' @param zoom_ylim Specific zoom ranges for y axis.
#' @param saveppt Output PPT format.
#' @param savehtml Save the images as HTML files.
#' @param ... Extra parameters to \code{\link[ggplot2]{ggsave}}.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' sp_ggplot_layout(p)
#'
#' ## End(Not run)
#'
sp_ggplot_layout <-
  function(p,
           xtics_angle = 0,
           legend.position = "right",
           extra_ggplot2_cmd = NULL,
           filename = NULL,
           x_label = NULL,
           y_label = NULL,
           title = NULL,
           coordinate_flip = FALSE,
           ylim = NULL,
           width=12,
           height=6.18,
           fontname = '',
           base_font_size = 10,
           additional_theme = NULL,
           zoom_split = FALSE,
           zoom_xlim = NULL,
           zoom_ylim = NULL,
           saveppt = FALSE,
           savehtml = FALSE,
           ...) {
    p <-
      p + theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        text = element_text(family = fontname, face = "plain",
                            colour = "black", size = base_font_size,
                            lineheight = 0.9,  hjust = 0.5,
                            vjust = 0.5, angle = 0,
                            margin = margin(), debug = FALSE),
        axis.line.x = element_line(
          size = 0.4,
          colour = "black",
          linetype = 'solid'
        ),
        axis.line.y = element_line(
          size = 0.4,
          colour = "black",
          linetype = 'solid'
        ),
        axis.ticks = element_line(size = 0.4,
                                  colour = "black")
      )

    if (xtics_angle != 0) {
      if (xtics_angle == 90) {
        p <- p + theme(axis.text.x =
                         element_text(
                           angle = xtics_angle,
                           hjust = 1,
                           vjust = 0.5
                         ))
      } else if (xtics_angle == 45) {
        p <- p + theme(axis.text.x =
                         element_text(
                           angle = xtics_angle,
                           hjust = 0.5,
                           vjust = 0.5
                         ))
      } else {
        p <- p + theme(axis.text.x =
                         element_text(
                           angle = xtics_angle,
                           hjust = 0.5,
                           vjust = 0.5
                         ))
      }
    }

    if (!sp.is.null(x_label)) {
      p <- p + xlab(x_label)
    }

    if (!sp.is.null(y_label)) {
      p <- p + ylab(y_label)
    }

    if (!sp.is.null(title)) {
      p <- p + labs(title = title)
    }

    p <- p + theme(legend.position = legend.position)

    #add additional ggplot2 supported commands

    if (!sp.is.null(extra_ggplot2_cmd)) {
      p <- p + eval(parse(text = extra_ggplot2_cmd))
    }

    if (coordinate_flip) {
      p <- p + coord_flip()
    }

    if (!sp.is.null(ylim)){
      p <- p + coord_cartesian(ylim = ylim)
    }


    additional_theme <- additional_theme[!sapply(additional_theme, sp.is.null)]
    if(length(additional_theme)>0){
      p <- p + do.call(theme, additional_theme)
    }

    # if (!sp.is.null(zoom_variable) && !sp.is.null(zoom_range)){
    #   p <- p + eval(parse(text = paste("facet_zoom(",zoom_axis,"=",zoom_variable,"==c(",zoom_range,"))")))
    # }
    if (!sp.is.null(zoom_xlim) || !sp.is.null(zoom_ylim)){
      p <- p + facet_zoom(xlim = zoom_xlim, ylim = zoom_ylim, split = zoom_split)
    }


    # output pictures
    if (sp.is.null(filename)) {
      return(p)
    } else{
      ggsave(p,
             filename = filename,
             units = c("cm"),
             width = width,
             height = height,
			 # added for abnormal pdf output
			 useDingbats = FALSE,
             ...)

	  cwd = getwd()
	  #print(cwd)
	  #print(filename)
	  #if(grepl("Cloud_Platform", cwd)){
	  #	cwd = "/var/www/html/Cloud_Platform//Cloud_Platform/public/"
	    # filename_ = basename(filename)
		#filename = paste0(cwd, filename)
		#print(filename)
	  #}
      if (saveppt){
	  # print(filename)
	  # print(dirname(filename))
	  # print(getSrcDirectory(function(x) {x}))
	  # print(dirname(sys.frame(1)$ofile))
	  # normalizePath(paste0(getwd(),dirname(filename),sep="/"))
      eoffice::topptx(p, filename = paste0(filename,".pptx"),
             width = width, height = height)
      }
      if (savehtml){
      plot_p <- plotly::ggplotly(p)
      htmlwidgets::saveWidget(as_widget(plot_p), paste0(filename,".index.html"))
      }
    }
    p
  }

#' Get the x, y limits of a ggplot2 plot
#'
#' @param p A ggplot2 object
#'
#' @return A list list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
#' @export
#'
#' @examples
#' ## Not run:
#' sp_get_ggplot_limits(p)
#'
#' ## End(Not run)
sp_get_ggplot_limits <- function(p) {
  # https://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object#
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax
  )
}

#' Return if unique values of two vectors are the same (order does not matter)
#'
#' @param x A vector
#' @param y A vector
#'
#' @return Logial value T or F
#' @export
#'
#' @examples
#'
#' value.identical(c('a','a','b','d'), c('d','d','a','b'))
#'
#' # TRUE
#'
value.identical <-
  function(x, y, treat_fully_contain_as_identical = F) {
    # numerical columns not used for compare
    if(is.numeric(x) || is.numeric(y)){
      return(FALSE)
    }
    x_unique = sort(unique(as.character(x)))
    y_unique = sort(unique(as.character(y)))
    all_ident = identical(x_unique, y_unique)
    if (all_ident) {
      return(all_ident)
    }
    if (treat_fully_contain_as_identical) {
      xy_intersect = intersect(x_unique, y_unique)
      return(identical(xy_intersect, x_unique) ||
               identical(xy_intersect, y_unique))
    } else {
      return(all_ident)
    }
  }

#' Detect pairs of columns with same unique values (order does not matter) in two dataframes.
#'
#' @param df1 Dataframe1
#' @param df2 Dataframe2
#' @param only_allow_one_match Default FALSE. This parameters is designed to get only one pair of matched columns
#' between two dataframes to supply as parameters for \link{merge} function (when TRUE).
#'
#' @return A dataframe containing names of matched columns. Or a vetor containing names of matched columns
#' when `only_allow_one_match` is `TRUE` and there do have one match.
#' @inheritParams value.identical
#' @export
#'
#' @examples
#'
#' vec1 <- data.frame(col1=c('a','a','b','d'), a=c(1,2,3,4))
#' vec2 <- data.frame(col2=c('d','d','a','b'), b=c(1,2,4,5),a=c(1,2,3,4))
#' get_matched_columns_based_on_value(vec1, vec2)
#'
#' #     match_1 match_2
#' # DF1    col1       a
#' # DF2    col2       a
#'
#' vec2 <- data.frame(col2=c('d','d','a','b'))
#' get_matched_columns_based_on_value(vec1, vec2)
#'
#' #     match_1
#' # DF1    col1
#' # DF2    col2
#'
#' get_matched_columns_based_on_value(vec1, vec2, only_allow_one_match = T)
#'
#' # "col1" "col2"
#'
#'
get_matched_columns_based_on_value <-
  function(df1,
           df2,
           only_allow_one_match = F,
           treat_fully_contain_as_identical = T) {
    if (length(df1) == 1) {
      df1['__extra_s_p_column__'] = '__extra_s_p_column__'
    }
    if (length(df2) == 1) {
      df2['__extra_s_p_column__'] = '__extra_s_p_column2__'
    }
    df1_rownames = rownames(df1)
    df2_rownames = rownames(df2)
    # print(sp.is.null(df1_rownames))
    # print(sp.is.null(df2_rownames))

    if(value.identical(df1_rownames, df2_rownames, treat_fully_contain_as_identical) &&
       # ignore default number row names
       !(value.identical(df2_rownames, as.character(1:length(df2_rownames))))){
      print(1273)
      return(c(0,0))
    }

    matches <-
      sapply(df2, function(x)
        sapply(
          df1,
          value.identical,
          x,
          treat_fully_contain_as_identical = treat_fully_contain_as_identical,
          simplify = T
        ),
        simplify = T)
    # print(matches)
    matches_index <-
      as.data.frame(which(matches == T, arr.ind = TRUE))
    # print(matches_index)
    df1_colnames <- colnames(df1)
    df2_colnames <- colnames(df2)
    matches_names <- apply(matches_index, 1,
                           function(x)
                             c(df1_colnames[x[1]], df2_colnames[x[2]]))
    matches_names <- as.data.frame(matches_names)
    # print(matches_names)
    if (nrow(matches_names) == 0) {
      stop(
        "No columns matched each other between given two data.frames. The program does not know which to return. Please check."
      )
    }
    matches_names_count <- length(matches_names)
    # print(matches_names_count)
    colnames(matches_names) <-
      paste0("match_", 1:matches_names_count)
    rownames(matches_names) <- c("DF1", "DF2")
    # print(matches_names)
    if (only_allow_one_match)
      if (matches_names_count > 1) {
        print(
          "Multiple pairs of columns matched each other between given two data.frames. The program does not know which to return. Please check."
        )
        return(as.vector(matches_names[, 1]))
      } else {
        # In case there are factors
        return(as.vector(matches_names[, 1]))
      }

    matches_names
  }

#' Check matched columns between two data frames and try to merge them.
#'
#' @param df1 Dataframe 1
#' @param df2 Dataframe 2
#' @param ... Extra parameters given to \code{\link[base]{merge}}.
#'
#' @return merged dataframe
#' @export
#'
#' @examples
#'
#' vec1 <- data.frame(col1=c('a','a','b','d'), a=c(1,2,3,6))
#' vec2 <- data.frame(col2=c('d','d','a','b'), b=c(1,2,4,5),a=c(1,2,3,4))
#' merge_data_with_auto_matched_column(vec1, vec2)
#'
merge_data_with_auto_matched_column <- function(df1, df2, ...){
  matched_column <-
    get_matched_columns_based_on_value(df1, df2,
                                       only_allow_one_match = T)

  data <-
    merge(df1, df2, by.x = matched_column[1], by.y = matched_column[2], suffixes = c("",".y"), ...)

  if(matched_column[2] != matched_column[1] || matched_column[1] != 0) {
    data[[matched_column[2]]] = data[[matched_column[1]]]
  }

  invisible(data)
}

