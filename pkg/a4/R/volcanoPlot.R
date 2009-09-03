setGeneric("volcanoPlot", function(x, y, pointLabels, ...){
      standardGeneric("volcanoPlot")
})


setMethod("volcanoPlot",
    signature(x = "tTest", 
        y = "missing",
        pointLabels = "missing"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL){
            
      logRatio <- x[,"logRatio"]
      pValue <- x[,"p"]
      pointLabels <- x[,"gSymbol"]
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub)
                  
})

setMethod("volcanoPlot",
    signature(x = "tTest", 
        y = "missing",
        pointLabels = "character"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL){
      
      logRatio <- x[,"logRatio"]
      pValue <- x[,"p"]
      
      if (length(pValue) != length(pointLabels))
        stop("'pointLabels' should have the same length as the number of rows of 'x'")
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub)
      
})


setMethod("volcanoPlot",
    signature(x = "limma", 
        y = "missing",
        pointLabels = "missing"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL){

      logRatio <- as.matrix(x@MArrayLM$coef)[, 2]
      pValue <- as.matrix(x@MArrayLM$p.value)[, 2]
      pointLabels <- x@geneSymbols
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, logTransformP = TRUE, # p values
          smoothScatter = smoothScatter, xlab = xlab, ylab = ylab, 
          main = main, sub = sub)

})

setMethod("volcanoPlot",
    signature(x = "limma", 
        y = "missing",
        pointLabels = "character"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10, smoothScatter = TRUE, 
        xlab = NULL, ylab = NULL, main = NULL, sub = NULL){
      
      logRatio <- as.matrix(x@MArrayLM$coef)[, 2]
      pValue <- as.matrix(x@MArrayLM$p.value)[, 2]
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, logTransformP = TRUE, # p values  
          smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub)      
})


setMethod("volcanoPlot",
    signature(x = "numeric", 
        y = "numeric",
        pointLabels = "character"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL){
      
      if ((length(x) != length(y)) | (length(x) != length(pointLabels)))
        stop("'x', 'y' and 'pointLabels' should have equal length")
      
      volcanoplotter(logRatio = x, pValue = y, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub)
})

setMethod("volcanoPlot",
    signature(x = "numeric", 
        y = "numeric",
        pointLabels = "missing"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL){
      
      if (length(x) != length(y))
        stop("'x' and 'y' should have equal length")
      
      if (is.null(names(x))){
        if (is.null(names(y))){
          stop(paste("nor 'x' nor 'y' have names that can be used to use as default 'pointLabels'\n",
                     "please make sure either 'x' or 'y' has names or, alternatively, \n",
                     "explicitly provide a 'pointLabels' argument", sep = "")) 
        } else {
          pointLabels <- names(y)
        } 
      } else {
        pointLabels <- names(x)
      }
      
      
      volcanoplotter(logRatio = x, pValue = y, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub)
    })


### workhorse function common to all volcanoPlot methods

volcanoplotter <- function(logRatio, pValue, pointLabels,
    topPValues = 10, topLogRatios = 10, logTransformP = TRUE,
    smoothScatter = TRUE, xlab = NULL, ylab = NULL, main = NULL, 
    sub = NULL, newpage = TRUE){
  ### checks                      
  if (!is.numeric(logRatio))
    stop("'logRatio' should be numeric")
  if (!is.numeric(pValue))
    stop("'pValue' should be numeric")
  # test if pValue is between 0 and 1
    if (any(pValue < 0 | pValue > 1))
      stop("'pValue' should be >= 0 and <= 1") # prevent from being already on log scale
  
  pVals <- if (logTransformP) -log10(pValue) else pValue 
  
  ### compute which points to label on the graph
  topLR <- order(abs(logRatio), decreasing = TRUE)[seq(length.out = topLogRatios)] 
  topP <- order(pVals, decreasing = TRUE)[seq(length.out = topPValues)]#logTransformP)[seq(length.out = topPValues)]
  pointsToLabel <- union(topP, topLR)
  
  ### set up graph
  if (newpage)
    grid.newpage()
  pvp <- plotViewport(c(5, 6, 5, 3))
  pushViewport(pvp)
  
  # compute maximum grobwidth
  
  tg <- textGrob(label = pointLabels[pointsToLabel], # TV: problematic sub _ratio 
      x = unit(logRatio[pointsToLabel], "native"),
      y = unit(pVals[pointsToLabel], "native"), 
      gp = gpar(cex = 0.65, col = "black"))
  
  maxLabelWidth <- max(grobWidth(tg))
  nMaxLabelWidth <- convertHeight(maxLabelWidth, "native", 
      valueOnly = TRUE)
  
  dvp <- dataViewport(xscale = range(logRatio) + c(-nMaxLabelWidth/2.2, nMaxLabelWidth/2.2),
      yscale = range(pVals))
  
  pushViewport(dvp)
  
  atPositionsY <- grid.pretty(current.viewport()$yscale)
  # atPositionsX <- grid.pretty(current.viewport()$xscale)
  xa <- xaxisGrob(name = "xa")# , at = atPositionsX)# , vp = pvp)
  
  grid.yaxis(name = "ya", at = atPositionsY, label = signif(10^(-atPositionsY), 1))
  
  # move the x axis down a bit
  moveUnit <- unit(-0.5, "char")
  xa <- editGrob(xa, edits = gEditList(
          gEdit("major", y = moveUnit),   # defaults to 0npc
          gEdit("ticks", y0 = moveUnit),  # defaults to 0npc
          gEdit("ticks", y1 = unit(-0.5, "lines") + moveUnit),   # defaults to -0.5lines
          gEdit("labels", y = unit(-1.5, "lines") + moveUnit)))  # defaults to -1.5lines 
  grid.draw(xa)
  
  # box(bty = "l", lwd = 1.5)
  dotColors <- if (smoothScatter){ 
        densCols(x = logRatio[-pointsToLabel], y = pVals[-pointsToLabel])
      } else {
        brewer.pal(9, "Blues")[4]
      } 
  grid.points(x = unit(logRatio[-pointsToLabel], "native"),
      y = unit(pVals[-pointsToLabel], "native"),
      pch = 20, gp = gpar(col = dotColors))
  
  grid.draw(tg)
  
  if (!is.null(main)){
    grid.text(label = main, y = unit(1, "npc") + unit(2, "lines"),
        gp = gpar(fontface = "bold"))
  }
  if (!is.null(sub)){
    grid.text(label = sub, y = unit(-4.25, "lines") + 0.5 * moveUnit)
  }
  if (!is.null(xlab)){
    grid.text(label = xlab, y = unit(-3, "lines") + 0.5 * moveUnit)
  }
  if (!is.null(ylab)){
    grid.text(label = ylab, x = unit(-4.5, "lines"), rot = 90)
  }          
}
