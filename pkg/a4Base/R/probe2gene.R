#' Auxiliary function for (currently) spectralMap
#' allowing the conversion of Affy probeset IDs to gene symbols
#' 
#' @author Tobias Verbeke
#' @param probesetIDS Affy probeset IDs
#' @param chipPkg string indicating the annotation package for the chip
probe2gene <- function(probesetIds, chipPkg){
  # DB-based annotation
  featureSymbols <- aafSymbol(probesetIds, chipPkg)
  return(getText(featureSymbols))
}
