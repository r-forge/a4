# simulate data for package testing and demonstration purposes
# nCols: number of samples
# nRows: number of features (genes)
# nEffectRows: number of differentially expressed features
# nNoEffectCols: number of samples for which the profile 
# of a differentially expressed feature will be
# set similar to the other class 


simulateData <- function(nCols = 40, nRows = 1000, nEffectRows = 5, nNoEffectCols = 5,
    betweenClassDifference = 1, withinClassSd = 0.5){
  ## response
  if (nCols %% 2 == 1) stop("'nCols' should be even")
  if (nNoEffectCols < 0) stop("'nNoEffectCols' should be positive (or zero)")
  # TODO: only works when even number of nCols !! TV (cf. infra as well)
  
  yData <- as.factor(rep(c("A", "B"), each = nCols / 2))
  
  ## no effect
  xData <- matrix(rnorm(nRows * nCols, mean = 8.73), ncol = nCols)
  colnames(xData) <- paste("Sample", 1:nCols, sep = "")
  rownames(xData) <- paste("Gene", 1:nRows, sep = ".") # change TV to prevent 
  # possible bug in new MLInterfaces version
  
  ## add effect to a certain number of rows
  xDataEffect <- matrix(0, nrow = nEffectRows, ncol = nCols)
  
  if (nNoEffectCols == 0){
    
    if (nEffectRows > 0) {
      
      ## create signal  
      for (i in seq(nEffectRows)){
        
        xDataEffect[i, ] <- rnorm(nCols, mean = 0, sd = withinClassSd) + 
            betweenClassDifference * as.numeric(yData)
      }
      xData[seq(nEffectRows), ] <- xDataEffect
      
    }  
  } else {
    ## create misbehaving samples
    
    yDataAnti <- factor(levels(yData)[(as.numeric(yData) %% 2) + 1], levels = levels(yData))
    yDataNoEffectCols <- yData
    yDataNoEffectCols[1:nNoEffectCols] <- yDataAnti[1:nNoEffectCols]
    
    if (nEffectRows > 0) {
      ## create signal  
      for (i in seq(nEffectRows)){
        
        xDataEffect[i, ] <- rnorm(nCols, mean = 0, sd = withinClassSd) + 
            betweenClassDifference * as.numeric(yDataNoEffectCols)
      }
      xData[seq(nEffectRows), ] <- xDataEffect
    }
  }  
  
  ### turn data into an ExpressionSet object
  info <- data.frame(yData, row.names = colnames(xData))
  colnames(info) <- "type"
  labels <- data.frame(type = "type")
  pd <- new("AnnotatedDataFrame", data = info, varMetadata = labels)
  # new("ExpressionSet", exprs = as.matrix(xData), phenoData = pd)
  eset.all <- new("ExpressionSet", exprs = as.matrix(xData), phenoData = pd) 
  # "ExpressionSet")  # will generate a lot of warnings
}




