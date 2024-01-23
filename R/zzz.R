.onAttach <- function(libname, pkgname){
  packageStartupMessage(paste0("Welcome to ImageGP package which is the base package for plot functions of https://www.bic.ac.cn/BIC. \n",
                               "This package does not require you install all depended packages, ",
                               "since one may not need all functions in this package.\n",
                               "However, when there is a message implying some functions are missing,",
                               "please install these packages manually.\n",
                               "Or more specially, check Plot.Rmd in vignettes first.\n")
  )
}
