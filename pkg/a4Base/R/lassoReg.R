# 
# Lasso using package 'glmnet'
#  both classification and regression are possible
#
###############################################################################

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




