setGeneric("histPvalue", function(object, ...){
      standardGeneric("histPvalue")
})    

setMethod("histPvalue", "limma",
    function(object, ...){
    
    nGenes <- length(object@geneSymbols)
    pValue <- topTable(object, n = nGenes)$P.Value
    
    histpvalueplotter(pValue = pValue)      
})

setMethod("histPvalue", "numeric",
    function(object, ...){
    histpvalueplotter(pValue = object, ...)      
})


histpvalueplotter <- function(pValue, addLegend = FALSE, ...){
  
  histOutput <- hist(pValue, 50, col = "skyblue", main = "")
  lengthHist <- length(histOutput$counts)
  meanNonDE <- mean(histOutput$counts[(lengthHist/2):lengthHist])
  abline(h = meanNonDE, col = 'goldenrod', lwd = 2)
  
  if (addLegend){
    legend("topright", bty = "n", 
        legend = paste(propDEgenes(pValue), "% DE genes", sep = ""),
        text.col = "goldenrod", cex = 1.2)
  }
}
