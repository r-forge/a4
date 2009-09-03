.onAttach <- function(libname, pkgname){
  options(error = NULL)
  message(paste("\na4 version ", packageDescription("a4")$Version, 
          "\n", sep = ""))
}
