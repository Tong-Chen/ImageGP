
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'
#   Generate DOC (Mac):              'Ctrl + Shift + Option + r'

sp_data_pass_normality_test <- function(data, group_variable=NULL, value=NULL){
  if(is.vector(data)){
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
shapiro.test2 <- function(data, threshold=1000){
  len_data <- length(data)
  if(len_data < 3){
    print("Too little data for normality test!")
    return(FALSE)
  } else if(len_data<=threshold){
    return(shapiro.test(data)$p.value>0.05)
  } else {
    print(paste0("For data with more than ", threshold, " samples, always return TRUE."))
    return(TRUE)
  }
}
