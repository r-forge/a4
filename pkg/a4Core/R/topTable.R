setGeneric("topTable", function(fit, n, ...){ # common to nlcv and (at least) a4Classif
      standardGeneric("topTable")
    })

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

