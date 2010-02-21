#
# Filtering on Intensity and Variance
# 
# Author: wtalloen
###############################################################################



#' Filtering on Intensity and Variance
#' @param object  
#' @param IntCutOff intensity cutoff...
#' @param IntPropSamples 
#' @param VarCutOff 
#' @returnType 
#' @return 
#' @export
filterVarInt <- function(object,
    IntCutOff = log2(100), # exclude low-expressed genes: genes with less than 6.6 on log2 scale 
    IntPropSamples = 0.25, # in more than 3/4 of the samples
    VarCutOff = 0.5){			 # exclude small-variance genes: genes with a InterQuartileRange smaller than 0.5
  f1 <- pOverA(IntPropSamples, IntCutOff)     
  f2 <- function(x) IQR(x) > VarCutOff   
  ff <- filterfun(f1, f2)
  selected <- genefilter(object, ff)
  esSel <- object[selected, ]
  return(esSel)
}
