setGeneric("histPvalue", function(object, ...){
      standardGeneric("histPvalue")
})    

setMethod("histPvalue", "limma",
    function(object, ...){
    
    nGenes <- length(object@geneSymbols)
    
    # currently default of coef = 2 is OK (as only limmaTwoLevels generates an object of class 'limma')
    pValue <- topTable(object, coef = 2, n = nGenes)$P.Value
    
    histpvalueplotter(pValue = pValue, ...)      
})

setMethod("histPvalue", "MArrayLM",
    function(object, coef, ...){
      
      if (missing(coef))
        stop("Please specify a 'coef' argument to select a coefficient for the topTable function used internally.")
      
      pValue <- topTable(object, coef = coef, n = nrow(object))$P.Value
      
      histpvalueplotter(pValue = pValue, ...)      
    })


setMethod("histPvalue", "numeric",
    function(object, ...){
    histpvalueplotter(pValue = object, ...)      
})


histpvalueplotter <- function(pValue, addLegend = FALSE, xlab = NULL, ylab = NULL, main = NULL, ...){
  
  mainTitle <- if (is.null(main)) "" else  main
  
  histOutput <- hist(pValue, 50, col = "skyblue", main = mainTitle, xlab = xlab, ylab = ylab, ...)
  lengthHist <- length(histOutput$counts)
  meanNonDE <- mean(histOutput$counts[(lengthHist/2):lengthHist])
  abline(h = meanNonDE, col = 'goldenrod', lwd = 2)
  
  if (addLegend){
    legend("topright", bty = "n", 
        legend = paste(propDEgenes(pValue), "% DE genes", sep = ""),
        text.col = "goldenrod", cex = 1.2)
  }
}
