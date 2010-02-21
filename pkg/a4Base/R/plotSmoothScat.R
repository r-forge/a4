# Genome wide plots with density coloring
# 
# Author: Willem Talloen, Suzy Van Sanden
###############################################################################



plotComb2Samples <- function(object, x, y,
		trsholdX = NULL, trsholdY = NULL,
		probe2gene = TRUE, ...){
	
	x <- exprs(object)[, as.character(x)]
	y <- exprs(object)[, as.character(y)]
	
	gSymbol <- featureData(object)$`SYMBOL`

	### determine points to label
	if (is.null(trsholdX) & is.null(trsholdY)){
		pointsToLabel <- NULL
	} else {
		XpointsToLabel <- if (!is.null(trsholdX)){
		  min(trsholdX) < x & x < max(trsholdX)
		} else {
			FALSE
		}
        YpointsToLabel <- if (!is.null(trsholdY)){
		  min(trsholdY) < y & y < max(trsholdY)
		} else {
			FALSE
		}
		pointsToLabel <- XpointsToLabel & YpointsToLabel
	}
	
	plot(x, y, axes=F, type="n", ...)
	axis(1, lwd = 1.5, las = 1); axis(2, lwd = 1.5, las = 1)
	box(bty='l',lwd = 1.5)
	
	if (!is.null(trsholdX) | !is.null(trsholdY)) {
		dotColors <- densCols(x[-pointsToLabel], y[-pointsToLabel])
		points(x[-pointsToLabel], y[-pointsToLabel],pch = 20, cex = 1, col = dotColors)
		text(x[pointsToLabel], y[pointsToLabel], labels = gSymbol[pointsToLabel],
				cex = 0.65, col = "black")
	} else {
		dotColors <- densCols(x, y)
		points(x, y,pch = 20, cex = 1, col = dotColors)
	}
}

#### Scatterplot matrix with density-dependent coloring
panel.plotSmoothScat <- function(x, y, ...) {
	points(x, y, axes=F, type="n", main="", xlab="", ylab="", ...)
	dotColors <- densCols(x, y)
	points(x, y, pch = 20, cex = 1, col = dotColors)
	abline(a=0, b=1, col="red")
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	
	test <- cor.test(x,y)
	# borrowed from printCoefmat
	Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
			cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
			symbols = c("***", "**", "*", ".", " "))
	
	text(0.5, 0.5, txt, cex = cex * r)
	text(.8, .8, Signif, cex=cex, col=2)
}

plotCombMultSamples <- function(exprsMatrix){
	pairs(exprsMatrix, lower.panel = panel.plotSmoothScat, upper.panel = panel.cor)
}
