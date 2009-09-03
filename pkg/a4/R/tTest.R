# 
# ordinary t-test
# 
# Author: wtalloen
###############################################################################

tTest <- function(object, groups, probe2gene = TRUE){

	# t-test for differential expression
	ttests <- rowttests(object, groups)
	pTtest <- data.frame(rownames(ttests), ttests[,'p.value'])
	statTtest <- data.frame(rownames(ttests), ttests[,'statistic'])
	
	# adjustment for multiple testing
	pAdjusted <- mt.rawp2adjp(ttests[, "p.value"], proc = c("BH"))
	pTtestBH <- pAdjusted$adjp[order(pAdjusted$index), "BH"]
	
	# log-ratio of differential expression
	labels <- as.numeric(factor(pData(object)[,groups]))-1
	logRatio <- rowMeans(exprs(object)[, labels == 1]) - rowMeans(exprs(object)[,labels == 0])
	
	if (probe2gene){
		gSymbol <- featureData(object)$`Gene Symbol`
		if (is.null(gSymbol))
			stop("There is no variable named'Gene Symbol' in the pData of the object.\n
				You may want to set the argument 'probe2gene' to FALSE (the default is TRUE)")
		
		pvalues <- data.frame(gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("gSymbol", "p", "logRatio", "pBH", "tStat")
	} else {
		pvalues <- data.frame(# gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("p", "logRatio", "pBH", "tStat")
	}
	
	class(pvalues) <- c("tTest", class(pvalues))
	return(pvalues)
}

fTest <- function(object, groups, probe2gene = TRUE, varEqual = FALSE){
	# t-test for differential expression
	ttests <- rowFtests(object, groups, var.equal = varEqual)
	pTtest <- data.frame(rownames(ttests), ttests[,'p.value'])
	statTtest <- data.frame(rownames(ttests), ttests[,'statistic'])
	
	# adjustment for multiple testing
	pAdjusted <- mt.rawp2adjp(ttests[, "p.value"], proc = c("BH"))
	pTtestBH <- pAdjusted$adjp[order(pAdjusted$index), "BH"]
	
	# log-ratio of differential expression
	labels <- as.numeric(factor(pData(object)[,groups]))-1
	logRatio <- rowMeans(exprs(object)[, labels == 1]) - rowMeans(exprs(object)[,labels == 0])
	
	if (probe2gene){
		gSymbol <- featureData(object)$`Gene Symbol`
		pvalues <- data.frame(gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("gSymbol", "p", "logRatio", "pBH", "fStat")
	} else {
		pvalues <- data.frame(# gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("p", "logRatio", "pBH", "fStat")
	}
	
	class(pvalues) <- c("fTest", class(pvalues))
	return(pvalues)
}

setOldClass("tTest")
setMethod("topTable", "tTest",
    function(fit, n){
      head(fit, n = n)
})

setOldClass("fTest")
setMethod("topTable", "fTest",
		function(fit, n){
			head(fit, n = n)
})

if (FALSE){ # TODO: remove
  require(plyr)
  x <- matrix(runif(970), ncol = 97)
  ttestfun <- function(x) t.test(x, var.equal = TRUE)
  aaply(x, ttestfun)
}

tTest2 <- function(object, groups, probe2gene = TRUE){
	groups <- pData(object)[, groups]
	testData <- exprs(object)
	# t-test for differential expression
	ttestfun <- function(y) t.test(y ~ groups, var.equal = TRUE)
	ttests <- aaply(testData, 1, ttestfun)
	rowttests(object, groups)
	pTtest <- data.frame(rownames(ttests), ttests[,'p.value'])
	statTtest <- data.frame(rownames(ttests), ttests[,'statistic'])
	
	# adjustment for multiple testing
	pAdjusted <- mt.rawp2adjp(ttests[, "p.value"], proc = c("BH"))
	pTtestBH <- pAdjusted$adjp[order(pAdjusted$index), "BH"]
	
	# log-ratio of differential expression
	labels <- as.numeric(factor(pData(object)[,groups]))-1
	logRatio <- rowMeans(exprs(object)[, labels == 1]) - rowMeans(exprs(object)[,labels == 0])
	
	if (probe2gene){
		gSymbol <- featureData(object)$`Gene Symbol`
		pvalues <- data.frame(gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("gSymbol", "p", "logRatio", "pBH", "tStat")
	} else {
		pvalues <- data.frame(# gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("p", "logRatio", "pBH", "tStat")
	}
	
	class(pvalues) <- c("tTest", class(pvalues))
	return(pvalues)
}
