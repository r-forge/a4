
#=======================
# plot1gene : make fixed or flexible (default) if requested

a4palette <- function(n){
  if (!is.numeric(n) | n < 1)
    stop("'n' should be a positive integer")
	res <- if (n==1) "red"
			else if (n==2) c("red","blue")
			else if (n==3) c("red","green","blue")
			else if (n==4) c("red","green","blue","purple")
			else rainbow(n)
  return(res)
}

### intensity plot
plot1gene <- function (probesetId = NULL, 
    geneSymbol = NULL, 
    object, groups, main = NULL, colvec = NULL,
    colgroups = NULL, 
    probe2gene = TRUE, sampleIDs = TRUE, 
    addLegend = TRUE, legendPos = "topleft", ...) {
  
  if (!is.null(probesetId) & !is.null(geneSymbol))
    stop("Please provide either a 'probesetId' or a 'geneSymbol'")
  
  if((length(sampleIDs) > 1) | !(is.logical(sampleIDs) | is.character(sampleIDs)))
    stop("'sampleIDs' should either be a logical or a character of length one")
  
  groups <- pData(object)[, groups]
  if (!is.factor(groups))
    stop("The variable referenced by 'groups' should be a factor")
  groups <- groups[, drop=TRUE] # remove unused groups
  nByGroup <- as.numeric(table(groups))
  
  if (is.null(geneSymbol)){ # probeset given
    probesetId <- as.character(probesetId)     # use names not position !!
    plotData <- exprs(object)[probesetId, ]
  } else { # gene given
    probesetPos <- which(geneSymbol == featureData(object)$`Gene Symbol`)
    if (!length(probesetPos))
      stop("gene 'gene' does not occur in ExpressionSet 'object'")
    
    probesetId <- featureNames(object)[probesetPos]
    if (length(probesetId) > 1)
      warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
              "probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
    plotData <- exprs(object)[probesetId[1], ] # use names not position !!
  }
  
  # order data by factor level
  orderGroups <- order(groups)
  groups <- groups[orderGroups]
  numericGroups <- as.numeric(groups) # after ordering
  
  plotData <- plotData[orderGroups]
  nc <- length(plotData)
  
  # define colors
  if (!is.null(colgroups)) {
    colgroups <- pData(object)[, colgroups]
    if (!is.factor(colgroups))
      stop("The variable referenced by 'colgroups' should be a factor")
    colgroups <- colgroups[, drop=TRUE] # remove unused groups
    colgroups <- colgroups[orderGroups]
    numericColgroups <- as.numeric(as.factor(colgroups)) # after ordering
    colGroupsNotDefinedAsArgument <- FALSE
  } else {
    colgroups <- groups
    numericColgroups <- numericGroups
    colGroupsNotDefinedAsArgument <- TRUE
  }
  
  if (is.null(colvec)){
    colvec <- a4palette(nlevels(colgroups))	  
  } else {
    if(length(colvec) != nlevels(colgroups))
      stop("'colvec' should contain as many elements as there are levels in 'groups' or 'colgroups'")
  }
  
  # prepare title
  if (probe2gene){
    gSymbol <- featureData(object)[probesetId[1],]$`Gene Symbol`
  }
  
  mainTitle <- if (is.null(main)){
        if (probe2gene)
          paste(gSymbol, " (", probesetId[1], ")", sep = "")
        else probesetId[1]
      } else { 
        main 
      }
  
  plot(1:(nc + 1), c(mean(plotData), plotData), type = "n",
      axes = FALSE, xlab = "", ylab = expression(log[2] ~ intensity),
      main = mainTitle)
  
  points(2:(nc + 1), plotData, bg = colvec[numericColgroups], pch = 21, cex=1.5)
  
  axis(2, las = 2, cex.axis = 0.7, lwd = 1.5)
  
  if (is.logical(sampleIDs)){
    if (sampleIDs){
      axis(1, las = 3, at = c(2:(nc + 1)), labels = names(plotData),
          cex.axis = 0.7, lwd = 1.5, las = 2)  
    } else {
      emptyLabels <- rep("", length(plotData))
      axis(1, las = 3, at = c(2:(nc + 1)), labels = emptyLabels,
          cex.axis = 0.7, lwd = 1.5, las = 2)
    }
  } else {
    sampleIDs <- pData(object)[, sampleIDs]
    axis(1, las = 3, at = c(2:(nc + 1)), labels = sampleIDs[orderGroups],
        cex.axis = 0.7, lwd = 1.5, las = 2)
  }
  mtext("labels Means", 1, at = 1, font = 2, cex = 0.7, las = 3)
  
  # add lines
  j <- 0
  if (colGroupsNotDefinedAsArgument) {
    for (i in 1:max(numericGroups)) {
      med <- mean(as.vector(plotData[(j + 1):(j + nByGroup[i])], mode = "numeric"))
      
      lines(x = c(0, 1.5), y = c(med, med), col = colvec[i], lty = 1, lwd = 2)
      lines(x = c(j + 2, j + 1 + nByGroup[i]), y = c(med, med), col = colvec[i], lty = 1, lwd = 2)
      
      j <- j + nByGroup[i] 
    }
  } else {			# if colored by an extra factor, the median lines should be in grey 
    for (i in 1:max(numericGroups)) {
      med <- mean(as.vector(plotData[(j + 1):(j + nByGroup[i])], mode = "numeric"))
      
      lines(x = c(0, 1.5), y = c(med, med), col = 'grey', lty = 1, lwd = 2)
      lines(x = c(j + 2, j + 1 + nByGroup[i]), y = c(med, med), col = 'grey', lty = 1, lwd = 2)
      
      j <- j + nByGroup[i] 
    } 
  } 
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = levels(colgroups),
        text.col = colvec, cex=1)
  }
  invisible(probesetId)
}

# plot boxplot per gene with raw data superimposed
boxPlot <- function(probesetId = NULL, 
    geneSymbol = NULL, 
    object, groups, main = NULL, colvec = NULL,
    colgroups = NULL, probe2gene = TRUE, 
    addLegend = TRUE, legendPos = "topleft", ...) {
  
  if (!is.null(probesetId) & !is.null(geneSymbol))
    stop("Please provide either a 'probeset' or a 'gene'")
  
  groups <- pData(object)[, groups]
  if (!is.factor(groups))
    stop("The variable referenced by 'groups' should be a factor")
  groups <- groups[, drop=TRUE] # remove unused groups
  numericGroups <- as.numeric(groups) # after ordering
  
  if (is.null(geneSymbol)){ # probeset given
    probesetId <- as.character(probesetId)     # use names not position !!
    plotData <- exprs(object)[probesetId, ]
  } else { # gene given
    probesetPos <- which(geneSymbol == featureData(object)$`Gene Symbol`)
    if (!length(probesetPos))
      stop("gene 'gene' does not occur in ExpressionSet 'object'")
    
    probesetId <- featureNames(object)[probesetPos] 
    if (length(probesetId) > 1)
      warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
              "probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
    plotData <- exprs(object)[probesetId[1], ] # use names not position !!
  }
  
  # define colors
  if (!is.null(colgroups)) {
    colgroups <- pData(object)[, colgroups]
    if (!is.factor(colgroups))
      stop("The variable referenced by 'colgroups' should be a factor")
    colgroups <- colgroups[, drop=TRUE] # remove unused groups
    numericColgroups <- as.numeric(as.factor(colgroups)) # after ordering
    colGroupsNotDefinedAsArgument <- FALSE
  } else {
    colgroups <- groups
    numericColgroups <- numericGroups
    colGroupsNotDefinedAsArgument <- TRUE
  }
  
  if (is.null(colvec)){
    colvec <- a4palette(nlevels(colgroups))	  
  } else {
    if(length(colvec) != nlevels(colgroups))
      stop("'colvec' should contain as many elements as there are levels in 'groups' or 'colgroups'")
  }
  
  # prepare title
  if (probe2gene){
    gSymbol <- featureData(object)[probesetId[1],]$`Gene Symbol`
  }
  
  mainTitle <- if (is.null(main)){
        if (probe2gene)
          paste(gSymbol, " (", probesetId[1], ")", sep = "")
        else probesetId
      } else { 
        main 
      }
  
  # make plot
  boxplot(plotData~groups, outline=FALSE, col='grey', type='n',
      ylab=expression(log[2]~concentration), las=1, main = mainTitle, ...)
  points(jitter(as.numeric(groups)),plotData,
      bg = colvec[numericColgroups], pch = 21, cex=1.5)
  
  # add legend
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = levels(colgroups),
        text.col = colvec, cex=1)
  }
  invisible(probesetId)
}

# plot combination of two genes
plotCombination2genes <- function(probesetId1=NULL, probesetId2=NULL,
    geneSymbol1=NULL, geneSymbol2=NULL,
    object, groups, addLegend = TRUE, legendPos = "topleft",
    probe2gene = TRUE, colvec = NULL, ...) {
  
  if (!is.null(probesetId1) & !is.null(probesetId2) & !is.null(geneSymbol1) & !is.null(geneSymbol2))
    stop("Please provide either two probesets or two genes")
  
  # change gene into probeset
  if (is.null(geneSymbol1)){ # probeset given
    probesetId1 <- as.character(probesetId1)     # use names not position !!
  } else { # gene given
    probesetPos1 <- which(geneSymbol1 == featureData(object)$`Gene Symbol`)
    if (!length(probesetPos1))
      stop("gene 'gene1' does not occur in ExpressionSet 'object'")
    
    probesetId1 <- featureNames(object)[probesetPos1]
    if (length(probesetId1) > 1)
      warning(paste("Gene1", geneSymbol1, "corresponds to", length(probesetId1), 
              "probesets; only the first probeset (",probesetId1[1],") has been displayed on the plot."))
  }
  
  if (is.null(geneSymbol2)){ # probeset given
    probesetId2 <- as.character(probesetId2)     # use names not position !!
  } else { # gene given
    probesetPos2 <- which(geneSymbol2 == featureData(object)$`Gene Symbol`)
    if (!length(probesetPos2))
      stop("gene 'gene2' does not occur in ExpressionSet 'object'")
    
    probesetId2 <- featureNames(object)[probesetPos2]
    if (length(probesetId2) > 1)
      warning(paste("Gene2", geneSymbol2, "corresponds to", length(probesetId2), 
              "probesets; only the first probeset (",probesetId2[1],") has been displayed on the plot."))
  }
  groups <- factor(pData(object)[, groups])[, drop=TRUE]
  exprGene1 <- exprs(object)[probesetId1[1], ]
  exprGene2 <- exprs(object)[probesetId2[1], ] 
  
  if (probe2gene){
    gSymbol1 <- featureData(object)[probesetId1[1],]$`Gene Symbol`
    gSymbol2 <- featureData(object)[probesetId2[1],]$`Gene Symbol`
  }
  
  if (is.null(colvec))
    colvec <- a4palette(nlevels(groups))
  
  # make plot
  plot(exprGene1, exprGene2, type = "n", 
      xlab = if (probe2gene) gSymbol1 else probesetId1[1], 
      ylab = if (probe2gene) gSymbol2 else probesetId2[1], 
      las = 1, ...)
  points(exprGene1, exprGene2, pch=21, cex=1.5,
      bg = colvec[as.numeric(groups)], ...)
  
  # add legend
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = levels(groups),
        text.col = colvec, cex=1)
  }
  invisible(list(probeset1=probesetId1, probeset2=probesetId2))
}

# parallel coordinate plots
profilesPlot <- function (object, probesetIds, sampleIDs = TRUE, 
    addLegend = TRUE, legendPos = "topleft", colvec = NULL,
    orderGroups = NULL, ...) {
  
  if (length(probesetIds) < 2)
    stop("Please provide at least two 'probesetIds'")
  
  if((length(sampleIDs) > 1) | !(is.logical(sampleIDs) | is.character(sampleIDs)))
    stop("'sampleIDs' should either be a logical or a character of length one")
  
  plotData <- t(exprs(object)[probesetIds,])
  
  # order data by factor level
  if (!is.null(orderGroups)){
    orderGroups <- as.factor(as.character(pData(object)[, orderGroups]))
    orderGroups2 <- order(orderGroups)
    plotData <- plotData[orderGroups2,]
  }
  
  if (is.null(colvec)){
    colvec <- a4palette(ncol(plotData))	  
  } else {
    if(length(colvec) != ncol(plotData))
      stop("'colvec' should contain as many elements as there are levels in 'groups' or 'colgroups'")
  }
  
  # plot
  matplot(plotData, type = "l",
      xlab = "", ylab = expression(log[2]~concentration),
      axes = FALSE, lwd = 1, lty = 1, col = colvec, ...)
  axis(2, las = 2, cex.axis = 0.7, lwd = 1.5)
  
  #x-axis
  if (is.logical(sampleIDs)){
    if (sampleIDs){
      axis(1, labels = rownames(plotData), las = 3, at = 1:nrow(plotData),
          cex.axis = 0.7, lwd = 1.5)
    } else {
      emptyLabels <- rep("", length(plotData))
      axis(1, labels = emptyLabels, las = 3, at = 1:nrow(plotData),
          cex.axis = 0.7, lwd = 1.5)
    }
  } else {
    sampleIDs <- pData(object)[, sampleIDs]
    axis(1, labels = sampleIDs[orderGroups2], las = 3, at = 1:nrow(plotData),
        cex.axis = 0.7, lwd = 1.5)
  }
  
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = colnames(plotData),
        text.col = colvec, cex=1)
  }
}
