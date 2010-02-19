confusionMatrix <- function(x, ...){
  UseMethod("confusionMatrix")
}

setGeneric("topTable", function(fit, n, ...){
      standardGeneric("topTable")
})

