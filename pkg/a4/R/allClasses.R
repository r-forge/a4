# TODO: Add comment
# 
# Author: tobias
###############################################################################

### define the annotationTable object
#' @author Tobias Verbeke
#' @param
setClass("annotationTable",
    representation = representation(displayData = "data.frame",
        displayCols = "list",
        hrefData = "data.frame"),
    prototype = list(displayData = data.frame(character()),
        displayCols = list(),
        hrefData = data.frame(character())))

### validity
.annotationTable.valid <- function(object){
  
  dimD <- dim(object@displayData)
  dimH <- dim(object@hrefData)
  
  if (!all.equal(dimD, dimH)) {
    warning("The displayData and hrefData should have the same dimensions")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

setValidity("annotationTable", .annotationTable.valid)


setClass("limma",
    representation = list(MArrayLM = "MArrayLM",
        geneSymbols = "character")
)

### S4 class for output of computeLogRatio function
setClass("ExpressionSetWithComputation",
		contains = "ExpressionSet")


