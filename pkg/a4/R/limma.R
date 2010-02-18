#' Wrapper for the limma function for the comparison of two groups
#' 
#' @param object object of class ExpressionSet
#' @param group string indicating the variable defining the two groups
#'               to be compared
limma2Groups <- function(object, group, probe2gene = TRUE){
  f <- factor(pData(object)[, group])[,drop = TRUE]
  if (nlevels(f) != 2)
    stop("Use 'limmaTwoGroups' only with a 'group' variable having two group levels")
  
  design <- model.matrix(~ f)
  fit <- lmFit(object, design)
  fit <- eBayes(fit)
  limmaObj <- new("limma", MArrayLM = fit, 
      geneSymbols = if (probe2gene) featureData(object)$`SYMBOL` else NULL)
      # use gene symbols from ExpressionSet and not the ones
      # that are in (the third column of) limmaObj@MArrayLM$genes
  
  return(limmaObj) 
}


setMethod("topTable", "limma",
    function(fit, n, coef = 2, genelist = fit$genes, 
        eb = fit[c("t", "p.value", "lods")], adjust.method = "BH",
        sort.by = "B", resort.by = NULL, p.value = 1, lfc = 0){
  
      # coef = 2 because we are not interested whether the intercept is significant 
      # but whether group 2 is significantly different from group 1
      
      fit <- fit@MArrayLM
      ### from limma:::topTable
      
      if (length(coef) > 1) {
        coef <- unique(coef)
        if (length(fit$coef[1, coef]) < ncol(fit)) 
          fit <- eBayes(fit[, coef])
        if (sort.by == "B") 
          sort.by <- "F"
        return(topTableF(fit, number = number, genelist = genelist, 
                adjust.method = adjust.method, sort.by = sort.by, 
                p.value = p.value))
      }
      fit <- unclass(fit)
      toptable(fit = fit[c("coefficients", "stdev.unscaled")], 
          coef = coef, number = n, genelist = fit$genes, A = fit$Amean, 
          eb = fit[c("t", "p.value", "lods")], adjust.method = "BH", 
          sort.by = sort.by, resort.by = resort.by, p.value = p.value, 
          lfc = lfc)
})

### redefine as well for MArrayLM objects

setMethod("topTable", "MArrayLM",
    function(fit, n, coef = 2, genelist = fit$genes, 
        eb = fit[c("t", "p.value", "lods")], adjust.method = "BH",
        sort.by = "B", resort.by = NULL, p.value = 1, lfc = 0){
      
      # coef = 2 because we are not interested whether the intercept is significant 
      # but whether group 2 is significantly different from group 1

      ### from limma:::topTable
      
      if (length(coef) > 1) {
        coef <- unique(coef)
        if (length(fit$coef[1, coef]) < ncol(fit)) 
          fit <- eBayes(fit[, coef])
        if (sort.by == "B") 
          sort.by <- "F"
        return(topTableF(fit, number = number, genelist = genelist, 
                adjust.method = adjust.method, sort.by = sort.by, 
                p.value = p.value))
      }
      fit <- unclass(fit)
      toptable(fit = fit[c("coefficients", "stdev.unscaled")], 
          coef = coef, number = n, genelist = fit$genes, A = fit$Amean, 
          eb = fit[c("t", "p.value", "lods")], adjust.method = "BH", 
          sort.by = sort.by, resort.by = resort.by, p.value = p.value, 
          lfc = lfc)
    })


