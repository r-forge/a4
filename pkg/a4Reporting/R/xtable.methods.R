xtable.topTableGlmnet <- function(x, caption = NULL, label = NULL, align = NULL, 
    digits = NULL, display = NULL, ...){
  
  xtable:::xtable.data.frame(x$topList, caption = caption, label = label, align = align,
      digits = digits, display = display, ...)
}

xtable.pamClassConfusionTable <- function(x, caption = NULL, label = NULL, align = NULL, 
    digits = NULL, display = NULL, ...){
  
  # TODO: some preprocessing
  xtable:::xtable.matrix(x, caption = caption, label = label, align = align,
      digits = digits, display = display, ...)
}

xtable.topTablePam <- function(x, ...){
  xtable(x$topList, ...)
}

xtable.topTableRfClass <- function(x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, ...){
  xtable(x$topList, ...)
}

