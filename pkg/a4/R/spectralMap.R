#' Spectral Map according to JnJ Standards
#' @author Tobias Verbeke
setGeneric("spectralMap", function(object, groups, ...){
      standardGeneric("spectralMap")
})    



#' Spectral Map according to JnJ Standards
#' @author Tobias Verbeke
#' @param object object of class \code{ExpressionSet}
#' @param groups string giving the name of the phenoData variable
#'               that defines the different groups to compare
#' @param makeLognormal boolean; if \code{TRUE}, the data will be backtransformed 
#'                   in order to obtain data of lognormal shape
#' @param plot.mpm.args list of arguments to pass to the \code{plot.mpm} function
#' @seealso \code{\link[mpm]{plot.mpm}}
#' @examples es <- simulateData()
#'           spectralMap(object = es, groups = "type", probe2gene = FALSE) 
setMethod("spectralMap",
  signature(object = "ExpressionSet", 
            groups = "character"),
  function(object, groups, makeLognormal = TRUE,
           mpm.args = list(row.weight = "mean", # mpmObject 
               col.weight = "constant", 
               logtrans = TRUE),
           plot.mpm.args = list(
             zoom = c(1,2),  # only these arguments are included that differ from plot.mpm defaults     
             label.tol = 10,  # please refer to ?plot.mpm for more information
             rot = c(-1, 1), 
             sub = "",
             lab.size = 0.85,
             col.group = pData(object)[, groups],
             # colors = c("orange1", "red", rainbow(length(unique(col.group)), start=2/6, end=4/6)),
             
             colors = c("wheat", # gene color (if no smoothScatter is used)
                        "black", # color for genes considered to be outlying 
                        a4palette(length(unique(pData(object)[, groups])))), # colors for the groups 
             col.size = 2,
             do.smoothScatter = TRUE),
          probe2gene = TRUE){
          expressionData <- exprs(object)
          chip <- annotation(object)
          chipAnnotationPkg <- paste(chip, "db", sep = ".")
          
          if (length(groups) > 1){
            stop("'groups' should be a string (character vector of length one)")
          }
          
          # plot.mpm.args$col.group <- pData(object)[, groups] # TV: no longer needed (added to default list) 
          mpmInput <- if (makeLognormal){ 
            data.frame(rownames(expressionData), 2^expressionData) 
          } else {
            data.frame(rownames(expressionData), expressionData)
          }
          mpmInput <- na.omit(mpmInput)
          mpm.args$data <- mpmInput
          plot.mpm.args$x <- do.call("mpm", mpm.args)    
          
          # adjust sample names (to escape the constraints of data frame column names)
          #   otherwise 'X' will have been prepended by the mpm function and displayed as such
          plot.mpm.args$x$col.names <- sampleNames(object)
          
          plot.mpm.args$zoom <- if (is.null(plot.mpm.args$zoom)) c(1,2) 
            else plot.mpm.args$zoom
          plot.mpm.args$label.tol <- if (is.null(plot.mpm.args$label.tol)) 10
            else plot.mpm.args$label.tol
          plot.mpm.args$rot <- if (is.null(plot.mpm.args$rot)) c(-1, 1)
            else plot.mpm.args$rot
          plot.mpm.args$sub <- if (is.null(plot.mpm.args$sub)) ""
            else plot.mpm.args$sub
          plot.mpm.args$lab.size <- if (is.null(plot.mpm.args$lab.size)) 0.85
            else plot.mpm.args$lab.size
          plot.mpm.args$col.group <- if (is.null(plot.mpm.args$col.group)) pData(object)[, groups]
            else plot.mpm.args$col.group
          plot.mpm.args$colors <- if (is.null(plot.mpm.args$colors)) 
              c("wheat", "black", rainbow(length(unique(pData(object)[, groups]))))
            else
              plot.mpm.args$colors
          plot.mpm.args$col.size <- if (is.null(plot.mpm.args$col.size)) 2
            else plot.mpm.args$col.size
          
          if (probe2gene){
            plot.mpm.args$labels <- pData(featureData(object))[plot.mpm.args$x$row.names,"Gene Symbol"]
			if (is.null(plot.mpm.args$labels))
				stop("There is no variable named'Gene Symbol' in the pData of the object.\n
								You may want to set the argument 'probe2gene' to FALSE (the default is TRUE)")
			
          }
          mpmPlot <- do.call("plot.mpm", plot.mpm.args)
          invisible(mpmPlot)
})
