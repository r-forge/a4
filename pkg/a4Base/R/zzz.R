.onAttach <- function(libname, pkgname){
  options(error = NULL)
  message(paste("\na4Base version ", packageDescription("a4Base")$Version, 
          "\n", sep = ""))
}
