# first stab at providing xtable methods for dataframes with hyperlinks
#                          for inclusion in Sweave reports
# Author: Tobias Verbeke, 2008-05-19
###############################################################################

# creation of objects of class annotationTable
# a 'double' data frame
### define the annotationTable object
#' @author Tobias Verbeke
#' @param
setClass("annotationTable",
    representation = representation(displayData = "data.frame",
        displayCols = "list",
        hrefData = "data.frame"),
    prototype = list(displayData = data.frame(character()),
        displayCols = list(),
        hrefData = data.frame(character())))

### validity
.annotationTable.valid <- function(object){
  
  dimD <- dim(object@displayData)
  dimH <- dim(object@hrefData)
  
  if (!all.equal(dimD, dimH)) {
    warning("The displayData and hrefData should have the same dimensions")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

setValidity("annotationTable", .annotationTable.valid)

##############################
### annotationTable object ###
##############################


# utility functions for automatic generation of hrefData
generateEntrezIdLinks <- function(x){
  # code from annaffy package
  url <- "http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=DetailsSearch&Term="
  if (!length(x))
    return(character(0))
  return(paste(url, x, sep = ""))
}

generateGOIdLinks <- function(x){
  # code from annaffy package
  url <- "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=query&query="
  
  if (!length(x))
    return(character(0))
  url <- paste(url, x, sep = "")
  for(i in 2:length(x))
    url <- paste(url, x, sep = "%0a")
  return(url)
}

# displayData: dataframe which should be displayed 
# hrefData: dataframe with hyperlinks corresponding to 
#           the elements of the displayData data frame
# displayCols: optional argument that can be passed to
#            automatically generate hrefData; displayCols
#            should be a list of named character vectors 
#            of length one such
#            as list(gid = "EntrezId", goid = "GOId"). Currently
#            only values EntrezId and GOId are supported

annotationTable <- function(displayData, 
    displayCols = NULL, 
    hrefData = NULL){

  if (!is.null(displayCols)){
    if (any(!(unlist(displayCols) %in% c("EntrezId", "GOId")))) 
      stop("'displayCols' should be a named list of strings 'EntrezId' or 'GOId'")
  }
  if (is.null(hrefData) & length(displayCols)){
    # create and populate hrefData
    hrefData <- as.data.frame(matrix("", nrow = nrow(displayData), 
            ncol = ncol(displayData)))
    names(hrefData) <- names(displayData)
    for (iCol in seq(along = displayCols)){
      
      iLinkCol <- which(names(displayData) %in% names(displayCols)[iCol])
      iLinkColType <- displayCols[[iCol]]
      iLinkColValues <- switch(iLinkColType, 
          EntrezId = generateEntrezIdLinks(displayData[, iLinkCol]),
          GOId = generateGOIdLinks(displayData[, iLinkCol]))
      hrefData[, iLinkCol] <- iLinkColValues                           
    }
  } else {
    origHrefData <- hrefData
    hrefData <- as.data.frame(matrix("", nrow = nrow(displayData), 
            ncol = ncol(displayData)))
    names(hrefData) <- names(displayData)
    if (any(!(names(origHrefData) %in% names(hrefData)))) 
      stop("The column names of 'hrefData' should be all present in the column names of 'displayData'")
    hrefData[, names(origHrefData)] <- origHrefData
  }
  # create new annotationTable object
  res <- new("annotationTable", displayData = displayData, 
      hrefData = hrefData)
  return(res)
}                


### define a particular show method
setMethod("show", "annotationTable", function(object){
      cat("annotationTable with hyperlink annotation for columns:\n")
      hlinkedCols <- which(!sapply(object@hrefData, function(x) all(is.na(x))))
      cat(paste(names(object@displayData)[hlinkedCols], collapse = ", "), "\n\n")
      print(object@displayData)
    })

### TODO as.data.frame [strip off all annotation information]                

#########################
### xtable formatting ###
#########################

# turn S3 generic to S4
setGeneric("xtable")

annotationTableSanitize <- function(str) {
  result <- str
  # only sanitize part that will be displayed i.e. content of second
  # curly braces \href{}{}
  isHref <- length(grep("^\\\\href\\{.*\\}\\{.*\\}$", result))
  if (isHref) displayPart <- gsub("(^\\\\href\\{(.*)\\}\\{).*\\}", "\\1", result)
  result <- gsub("^\\\\href\\{(.*)\\}\\{(.*)\\}", "\\2", result)
  result <- gsub("\\\\","SANITIZE.BACKSLASH",result)
  result <- gsub("$","\\$",result,fixed=TRUE)
  result <- gsub(">","$>$",result,fixed=TRUE)
  result <- gsub("<","$<$",result,fixed=TRUE)
  result <- gsub("|","$|$",result,fixed=TRUE)
  result <- gsub("{","\\{",result,fixed=TRUE)
  result <- gsub("}","\\}",result,fixed=TRUE)
  result <- gsub("%","\\%",result,fixed=TRUE)
  result <- gsub("&","\\&",result,fixed=TRUE)
  result <- gsub("_","\\_",result,fixed=TRUE)
  result <- gsub("#","\\#",result,fixed=TRUE)
  result <- gsub("^","\\verb|^|",result,fixed=TRUE)
  result <- gsub("~","\\~{}",result,fixed=TRUE)
  result <- gsub("SANITIZE.BACKSLASH","\\",result,fixed=TRUE)
  if (isHref) result <- paste(displayPart, result, "}", sep = "")
  return(result)
}

# TODO: add methods for other signatures

setMethod("xtable",
  signature(x="annotationTable", caption = "missing", label = "missing", align = "missing",
      digits = "missing", display = "missing"), # for an object of class ...
  function(x, caption, label, align, digits, display){
    # strategy: construct new data frame with hyperlinks
    # and use classical xtable on this one
    d <- x@displayData
    h <- x@hrefData
    at <- character(prod(dim(d))) # annotationTable
    dim(at) <- dim(d)
    colnames(at) <- names(d)
    numericCols <- which(sapply(d, is.numeric))
    
    for (iCol in seq(ncol(d))){
      
      d[[iCol]] <- as.character(d[[iCol]])
      h[[iCol]] <- as.character(h[[iCol]])
      
      if (all(h[[iCol]] == "") | all(is.na(h[[iCol]]))){ # no hyperlinks
        at[, iCol] <- d[[iCol]]    # display data only, not transformed into character
      } else {
        tm <- matrix(c(h[[iCol]], d[[iCol]]), ncol = 2)
        at[, iCol] <- paste("\\", "href{", 
            apply(tm, 1, paste, collapse = "}{"), "}", sep = "")
      }
    }
    at <- as.data.frame(at, stringsAsFactors = FALSE)
    for (iCol in numericCols)
      at[, iCol] <- as.numeric(as.character(at[, iCol])) # quick and dirty
    # \href{http://www.wikibooks.org}{wikibooks home}
    xt <- xtable(at)
    class(xt) <- c("xtableAnnotationTable", class(xt))
    xt
})

setMethod("xtable",
    signature(x="annotationTable", caption="ANY", label="ANY", 
        align="ANY", digits = "ANY", display = "ANY"), # for an object of class ...
    function(x, caption, label, align, digits, display){
      # strategy: construct new data frame with hyperlinks
      # and use classical xtable on this one
      d <- x@displayData
      h <- x@hrefData
      at <- character(prod(dim(d)))
      dim(at) <- dim(d)
      colnames(at) <- names(d)
      numericCols <- which(sapply(d, is.numeric))
      
      for (iCol in seq(ncol(d))){
        
        d[[iCol]] <- as.character(d[[iCol]])
        h[[iCol]] <- as.character(h[[iCol]])
        
        if (all(h[[iCol]] == "") | all(is.na(h[[iCol]]))){ # no hyperlinks
          at[, iCol] <- d[[iCol]]    # display data only, not transformed into character
        } else {
          tm <- matrix(c(h[[iCol]], d[[iCol]]), ncol = 2)
          at[, iCol] <- paste("\\", "href{", 
              apply(tm, 1, paste, collapse = "}{"), "}", sep = "")
        }
      }
      at <- as.data.frame(at, stringsAsFactors = FALSE)
      for (iCol in numericCols)
        at[, iCol] <- as.numeric(as.character(at[, iCol])) # quick and dirty
      # \href{http://www.wikibooks.org}{wikibooks home}
      xt <- xtable(at, caption = caption,
          label = label, align = align, digits = digits, display = display)
      class(xt) <- c("xtableAnnotationTable", class(xt))
      xt
})

setMethod("xtable",
    signature(x="annotationTable", caption="ANY", label="ANY", 
        align="ANY", digits = "numeric", display = "ANY"), # for an object of class ...
    function(x, caption, label, align, digits, display){
      # strategy: construct new data frame with hyperlinks
      # and use classical xtable on this one
      d <- x@displayData
      h <- x@hrefData
      at <- character(prod(dim(d)))
      dim(at) <- dim(d)
      colnames(at) <- names(d)
      numericCols <- which(sapply(d, is.numeric))
      
      for (iCol in seq(ncol(d))){
        
        d[[iCol]] <- as.character(d[[iCol]])
        h[[iCol]] <- as.character(h[[iCol]])
        
        if (all(h[[iCol]] == "") | all(is.na(h[[iCol]]))){ # no hyperlinks
          at[, iCol] <- d[[iCol]]    # display data only, not transformed into character
        } else {
          tm <- matrix(c(h[[iCol]], d[[iCol]]), ncol = 2)
          at[, iCol] <- paste("\\", "href{", 
              apply(tm, 1, paste, collapse = "}{"), "}", sep = "")
        }
      }
      at <- as.data.frame(at, stringsAsFactors = FALSE)
      for (iCol in numericCols)
        at[, iCol] <- as.numeric(as.character(at[, iCol])) # quick and dirty
      # \href{http://www.wikibooks.org}{wikibooks home}
      xt <- xtable(at, caption = caption,
          label = label, align = align, digits = digits, display = display)
      class(xt) <- c("xtableAnnotationTable", class(xt))
      xt
    })

### print method    
print.xtableAnnotationTable <- function(x, ...){
   # xtable:::print.xtable(x, sanitize.text = annotationTableSanitize, ...)
   NextMethod(x, sanitize.text = annotationTableSanitize, ...)
} 

    