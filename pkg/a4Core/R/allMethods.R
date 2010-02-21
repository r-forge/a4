confusionMatrix <- function(x, ...){ # common to a4Classif (pamClass) and nlcv
  UseMethod("confusionMatrix")
}

setGeneric("topTable", function(fit, n, ...){ # common to nlcv and (at least) a4Classif
  standardGeneric("topTable")
})


