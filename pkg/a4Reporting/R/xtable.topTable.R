xtable.topTableGlmnet <- function(x, caption = NULL, label = NULL, align = NULL, 
    digits = NULL, display = NULL, ...){
  
  xtable:::xtable.data.frame(x$topList, caption = caption, label = label, align = align,
      digits = digits, display = display, ...)
}

