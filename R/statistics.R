


# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'
#   Generate DOC (Mac):              'Ctrl + Shift + Option + r'

sp_data_pass_normality_test <-
  function(data,
           group_variable = NULL,
           value = NULL) {
    if (is.vector(data)) {
      shapiro.test2(data)
    }
    # library(rstatix)
    # data %>% group_by(combine__grp__for__statistis_sp) %>% shapiro_test(vars="Expr")
  }

#' An inner version of `shapiro.test` with two modifications.
#'
#' * Always return FALSE for sample size less than 3.
#' * Always return TRUE for sample size larger than given threshold (default 1000)
#'
#' @param data A vector of observation values.
#' @param threshold A number to define **large** sample size. Default `1000`. For data with
#' more than given `threshold` samples, always assume they past normality test,
#'
#' @return A logical value
#' @export
#'
#' @examples
#'
#' shapiro.test2(c(1,2,3,4))
#'
shapiro.test2 <- function(data, threshold = 1000) {
  len_data <- length(data)
  if (len_data < 3) {
    print("Too little data for normality test!")
    return(FALSE)
  } else if (len_data <= threshold) {
    return(shapiro.test(data)$p.value > 0.05)
  } else {
    print(paste0(
      "For data with more than ",
      threshold,
      " samples, always return TRUE."
    ))
    return(TRUE)
  }
}

#' Perform statistical test using given methods.
#'
#' @param data A data matrix
#' @param stat_value_variable The column represents the statistical value information.
#' @param stat_group_variable The column represents the statistical group information.
#' @param statistical_method Statistical method. For two groups, default <t.test>. For more than two groups, default <aov>.
#' @param add_y Add positions for each statistical label.
#' @param statistical_threshold_for_letters Threshold for treating as significance, default 0.05.
#' @return A list. list(data=data, a data frame with statistical information, Tukey_HSD=Tukey_HSD)
#' @export
#'
#' @examples
#'
#‘ data <-
#'
sp_diff_test <-
  function(data,
           stat_value_variable,
           stat_group_variable,
           statistical_method = "aov",
           add_y = TRUE,
           statistical_threshold_for_letters = 0.05) {
    # 代码修改自 amplicon包 microbiota/amplicon
    # https://github.com/microbiota/amplicon/blob/master/R/alpha_boxplot.R
    data <- droplevels(data)

    formula = as.formula(paste(stat_value_variable, "~", stat_group_variable))
    if (statistical_method == "aov") {
      model = aov(formula, data = data)
      Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)

      Tukey_HSD_table = as.data.frame(Tukey_HSD[[stat_group_variable]])
      # print(Tukey_HSD_table)

      if (length(unique(data[[stat_group_variable]])) == 2) {
        library(agricolae)
        out = LSD.test(model, stat_group_variable, p.adj = "none", alpha=statistical_threshold_for_letters)
        #print(out)
        LSD.test_table = as.data.frame(out$statistics)
        stat = out$groups
        data$stat = stat[as.character(data[[stat_group_variable]]),]$groups
      } else{
        if (length(unique(data[[stat_group_variable]])) == 2) {
          Tukey.levels = Tukey_HSD[[stat_group_variable]][, 4, drop = F]
        } else{
          Tukey.levels = Tukey_HSD[[stat_group_variable]][, 4]
        }

        # print(Tukey.levels)
        Tukey.labels = data.frame(multcompView::multcompLetters(
                      Tukey.levels, threshold=statistical_threshold_for_letters)['Letters'])
        Tukey.labels$group = rownames(Tukey.labels)
        Tukey.labels = Tukey.labels[order(Tukey.labels$group),]
        data$stat = Tukey.labels[as.character(data[[stat_group_variable]]), ]$Letters
      }
    }


    max_y = sapply(split(data, data[[stat_group_variable]]),
                     function (x) {max(x[[stat_value_variable]])})
    #print(max_y)

    data$y = max_y[data[[stat_group_variable]]] * 1.06

    return(list(data=data, Tukey_HSD_table=Tukey_HSD_table))
  }

#' Perform statistical test using given methods for split data.
#'
#' @param group_variable The column represents the group information.
#' The data would be split by this group and diff test would be performed within each group.
#' @inheritParams sp_diff_test
#' @param ... Other paramters givene to \link[sp_diff_test].
#'
#' @return A list. list(data=data, a data frame with statistical information, Tukey_HSD=Tukey_HSD)
#' @export
#'
#' @examples
sp_multiple_group_diff_test <-
  function(data,
           stat_value_variable,
           stat_group_variable,
           group_variable = NULL,
           ...) {

    if(sp.is.null(group_variable)){
      dataList = sp_diff_test(data, stat_value_variable=stat_value_variable,
                                         stat_group_variable=stat_group_variable,
                                         ...)
    } else {
      dataTmpList <- lapply(split(data, data[[group_variable]]),
                            sp_diff_test,
             stat_value_variable=stat_value_variable,
             stat_group_variable=stat_group_variable,
             ...)

      data = do.call(rbind, lapply(dataTmpList, function (x) x$data))
      Tukey_HSD_table = do.call(rbind, lapply(dataTmpList, function (x) x$Tukey_HSD_table))
      dataList = list(data=data, Tukey_HSD_table=Tukey_HSD_table)
    }
  return(dataList)
  }


# print(data)
