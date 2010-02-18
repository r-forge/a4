# TODO: Add comment
# 
# Author: wtalloen
###############################################################################


logReg <- function(object, groups, probesetId = NULL, 
		geneSymbol = NULL, main = NULL, probe2gene = TRUE, ...){
	
	dotCol <- 'goldenrod'
	lineCol <- 'blue'
	
	if (!is.null(probesetId) & !is.null(geneSymbol))
		stop("Please provide either a 'probeset' or a 'gene'")
	
	### create gene expression vector
	if (is.null(geneSymbol)){ # probeset given
		probesetId <- as.character(probesetId)     # use names not position !!
		exprGene <- exprs(object)[probesetId, ]
	} else { # gene given
		probesetPos <- which(geneSymbol == featureData(object)$`SYMBOL`)
		if (!length(probesetPos))
			stop("gene 'gene' does not occur in ExpressionSet 'object'")
		
		probesetId <- featureNames(object)[probesetPos]
		if (length(probesetId) > 1)
			warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
							"probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
		exprGene <- exprs(object)[probesetId[1], ] # use names not position !!
	}
	
	labels <- as.factor(as.character(pData(object)[,groups]))
	
	# prepare title
	if (probe2gene){
		gSymbol <- featureData(object)[probesetId[1],]$`SYMBOL`
	}
	
	mainTitle <- if (is.null(main)){
				if (probe2gene)
					paste(gSymbol, " (", probesetId[1], ")", sep = "")
				else probesetId[1]
			} else { 
				main 
			}
	
	logResGene <- glm(labels ~ exprGene, family = 'binomial')
	logRes <- data.frame(x = exprGene, y = labels, fit = fitted(logResGene))
	logRes$y <- labels[]
	logRes <- logRes[order(logRes$x),]
	
	par(mar = c(5, 4, 4, 7) + 0.1)
	plot(logRes$x, as.numeric(logRes$y)-1, axes=FALSE,
			pch = 21, bg = dotCol, 
			xlab=expression(log[2] ~ intensity), 
			ylab= paste('Probability of being',levels(labels)[2]),
			main = mainTitle, ...)
	lines(logRes$x, logRes$fit, col = lineCol, lwd=2)
	axis(1, las=1); axis(2, las=1); box(bty='l') 
	axis(4, at= c(0,1), labels= levels(labels), 
			col.axis='black', tick = FALSE, las=1)
	par(mar = c(5, 4, 4, 2) + 0.1)
	
	return(logRes)
	#-----------
}

probabilitiesPlot <- function(proportions,
		classVar,
		sampleNames,
		plot = TRUE,
		barPlot = FALSE, # barplot or MCREstimate-type scores plot
		layout = TRUE, 
		main = NULL,
		sub = NULL,
		...){ # additional arguments such as main, sub etc.
	
	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	
	plotData <- proportions # vector of proportions
	names(plotData) <- sampleNames
	classVar <- as.factor(classVar) # after using names
	classVarLevels <- levels(classVar)
	
	##  data values are grouped according to their observed labels
	##  (e.g. responders / non responders)
	plotData <- plotData[order(classVar)] 
	classVar <- classVar[order(classVar)] # order the labels as well
	
	if (plot){
		if (!barPlot) {
			### layout
			if (layout) layout(matrix(1:2, ncol = 1), heights = c(6, 1))
			
			### upper plot
			plot(x = seq(along = plotData), # indices (names) of the samples 
					y = plotData,   # vector of proportions of misclassification for each sample
					ylim = c(0, 1),
					type = "n", las =  3, ann = FALSE, axes = FALSE, ...)
			
			axis(1, las = 3, at = seq(along = plotData), 
					labels = names(plotData), cex.axis = 0.5)
			axis(2, las = 2)
			
			# draw grid first...
			abline(h = 0.5, col = "grey")
			abline(v = seq(length(plotData)), lty = "dashed", col = "grey")
			
			# ... then add dots  
			points(x = seq(along = plotData), y = plotData,
					pch = ifelse(plotData >= 0.5, 19, 17), cex = 1,
					col = ifelse(plotData >= 0.5, "darkblue", "orange"))
			
			title(main = if(is.null(main)){ 
								"Proportions Plot"
							} else {
								main
							})
			
			### add observed class membership
			nClasses <- length(classVarLevels)
			classColors <- brewer.pal(nClasses+3,"YlGn")[2:(nClasses+1)] # from MCREstimate
			
			sampleLocations <- seq(along = classVar)
			rect(xleft = sampleLocations-0.5, ybottom = -0.5, 
					xright = sampleLocations + 0.5, ytop = -0.015, 
					col = classColors[as.numeric(classVar)], 
					border = classColors[as.numeric(classVar)])
			abline(v = sampleLocations, lty = 2, col = "grey")
			
			### lower plot (with legends)
			op <- par(mar = c(0,4,0,2))
			plot(c(0, 1), type = "n", ann = FALSE, axes = FALSE)
			legend("left", legend = classVarLevels, fill = classColors,
					bty = "n")
			legend("right", legend = c("0.5 <= score <=   1", 
							"   0 <= score <  0.5"), 
					pch = c(19, 17), pt.cex = 1.5, col = c("darkblue", "orange"), 
					bty = "n")
			par(op)
						
		} else {
			barplot(height = plotData, col = ifelse(plotData >= 0.5, "green", "red"),  
					las = 3, ...)
			abline(h = 0.5, col = "grey")    
		}
	}
	
	par(def.par)#- reset to default
	invisible(plotData) # named vector of proportion correctly classified
}

