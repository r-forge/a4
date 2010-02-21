### utility function to generate an ExpressionSet object
### from expression data and phenodata
### 

`createExpressionSet` <-
		function(
				exprs = new("matrix"),
				phenoData = new("AnnotatedDataFrame"),
				varMetadata= NULL,
				dimLabels = c("rowNames","colNames"),
				featureData = NULL, 
				experimentData = new("MIAME"), 
				annotation = character(0), 
				changeColumnsNames = TRUE,...){  
	

	if (nrow(phenoData) != ncol(exprs)) 
    stop('phenoData must have the same number of rows than exprs number of columns')  

	if (all(rownames(phenoData) != colnames(exprs))){
		stop("rownames of phenoData are not identical to colnames of exprs.\n")
	}
	
	if (!is(phenoData, "AnnotatedDataFrame")){
		# we must prepare AnnotatedDataFrame
		# check varMetadata consistency with phenoData matrix
		if (!is.null(varMetadata)){
			if ((nrow(varMetadata)!=ncol(phenoData)) | (!'labelDescription' %in% colnames(varMetadata))){
				warning(paste("varMetadata not compliant with phenoData", "maybe there is not a column called 'labelDescription'","check ?AnnotatedDataFrame","--- we will not use it",sep='\n'))
				phenoData <- new('AnnotatedDataFrame',
						data=phenoData,dimLabels=dimLabels)
			}
			else {
				phenoData <-new('AnnotatedDataFrame',
						data=phenoData,
						varMetadata=varMetadata,
						dimLabels=dimLabels)
			}
		} else phenoData <- new('AnnotatedDataFrame',
					data=phenoData,
					dimLabels=dimLabels)
	}
	
	if (changeColumnsNames){
		oldcolnames <- colnames(exprs)
		if ("colNames" %in% colnames(phenoData)){
			if (any(duplicated(phenoData$colNames))) {
				warnings("Cant' use 'colNames' as new colnames as it has duplicates -- we use V1-Vn new names")
				colnames(exprs) <- paste('V',1:ncol(exprs),sep='')
			}
			else 
				colnames(exprs) <- phenoData$colNames
		}
		if ((!"colNames" %in% colnames(phenoData))){
			newcolnames <- do.call('paste',c(as.list(pData(phenoData)),sep='.'))
			if (any(duplicated(newcolnames))) {
				newcolnames <- paste(newcolnames,replicates(newcolnames),sep='.')            
				colnames(exprs) <- newcolnames
			}
		}
		
		phenoData$.oldcolnames <- oldcolnames
	}
	rownames(pData(phenoData)) <- colnames(exprs)  
	if (is(exprs,'data.frame')) exprs <- as.matrix(exprs)
	if (!is.null(featureData)){
		out <- new('ExpressionSet', 
				exprs=exprs,
				phenoData = phenoData,
				featureData=featureData, 
				experimentData = experimentData, 
				annotation = annotation)
	}
	else {
		out <- new('ExpressionSet', 
				exprs=exprs,
				phenoData = phenoData,
				experimentData = experimentData, 
				annotation = annotation)
	}
	return(out)
}


#createExprSet <- function(ExprData, PhenoData){
#	require(affy)
#	if(nrow(PhenoData) != ncol(ExprData)){
#		stop("The number of rows in PhenoData is not equal to the number of 
#			columns in the ExprData.\n")
#	}
#	
#	if(all(rownames(PhenoData) != colnames(ExprData))){
#		stop("rownames of PhenoData are not identical to colnames of ExprData.\n")
#	}
#	myExprSet <- new("ExpressionSet", exprs = ExprData)
#	pData(myExprSet) <- PhenoData
#	return(myExprSet)
#}


### utility function to combine two ExpressionSet objects
### 

`combineTwoExpressionSet` <-
		function(x,y){
# prioritary  keep information from x, append assayData and phylo data from y
	out <- x
	outAssayData <- new.env()
	assign("exprs",cbind(assayData(x)$exprs,assayData(y)$exprs),
			envir=outAssayData)
	assayData(out) <- outAssayData
	outPhenoData <- new("AnnotatedDataFrame",
			data=rbind(pData(x),pData(y)),
			varMetadata = varMetadata(x)) 
	phenoData(out) <- outPhenoData
	return(out)
}
