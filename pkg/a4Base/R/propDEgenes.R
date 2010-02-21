setGeneric("propDEgenes", function(object, ...){
      standardGeneric("propDEgenes")
    })

setMethod("propDEgenes", "limma",
    function(object, ...){
      
      nGenes <- length(object@geneSymbols)
      pValue <- topTable(object, n = nGenes)$P.Value
      
      propdegenescalculation(pValue = pValue)
})

setMethod("propDEgenes", "numeric",
    function(object, ...){
      
      propdegenescalculation(pValue = object)
      
    })

propdegenescalculation <- function(pValue){
  NbDEgenes <- length(pValue) - (sum(pValue > 0.5)*2)
  return( round((100/ length(pValue) )* NbDEgenes, 1) )
}

