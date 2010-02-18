# 
# Lasso using package 'glmnet'
#  both classification and regression are possible
#
###############################################################################

lassoClass <- function(object, groups){
  labels <- factor(pData(object)[,groups])
  object <- object[,!is.na(labels)]
  labels <- factor(pData(object)[,groups])
  
  fit <- glmnet(t(exprs(object)), labels, family="binomial", alpha = 1)
  fit$featureData <- pData(featureData(object))
  return(fit)	
}

lassoReg <- function(object, covariate){
  covariateVector <- pData(object)[, covariate]
  if (!is.numeric(covariateVector))
	  stop("The argument 'covariate' needs to refer to CONTINUOUS variable
		from the 'phenoData(object)'")

  object <- object[,!is.na(covariateVector)]
  covariateVector <- pData(object)[,covariate]

  fit <- glmnet(t(exprs(object)), covariateVector, family="gaussian", alpha = 1)
  fit$featureData <- pData(featureData(object))
  return(fit) 
}

setOldClass("glmnet")

setMethod("topTable",
    "glmnet",
    function(fit, n){
      summary.output <- summary(fit)
      coef.output <- coef(fit) # extract coefficients at a single value of lambda
      last.coef.output <- coef.output[, ncol(coef.output)]
      selProbeSets <- last.coef.output[which(last.coef.output != 0)]
      selProbeSetsGeneSymbol <- fit$featureData[names(selProbeSets), "SYMBOL"]
      
      selGenesOutput <- cbind.data.frame(selProbeSetsGeneSymbol, selProbeSets)
      rownames(selGenesOutput) <- names(selProbeSets)
      colnames(selGenesOutput) <- c('Gene','Coefficient')
      # remove the estimate of the intercept (typically, but not always, the first row)
      exclIntercept <- which(rownames(selGenesOutput)%in%c('','(Intercept)'))
      selGenesOutput <- selGenesOutput[-exclIntercept,]
      
      numberSelGenes <- nrow(selGenesOutput)
      topList <- selGenesOutput[order(abs(selGenesOutput[,2]),decreasing=TRUE),][1:(min(n, numberSelGenes)),] # first row is the estimate of the intercept.
      res <- list(topList = topList, numberSelGenes = numberSelGenes, n = n)
      class(res) <- "topTableGlmnet"
      return(res)
    }
)

print.topTableGlmnet <- function(x,  ...){
  cat("The lasso selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n", sep = "")
  print(x$topList, ...)
}

xtable.topTableGlmnet <- function(x, caption = NULL, label = NULL, align = NULL, 
    digits = NULL, display = NULL, ...){
  
  xtable:::xtable.data.frame(x$topList, caption = caption, label = label, align = align,
      digits = digits, display = display, ...)
}

