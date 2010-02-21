

### utility function to transform an ExpressionSet into an ExpressionSet with meta data similar
### to the ExpressionSets used in the pipeline

addGeneInfo <- function(eset, annotationLibrary = NULL){
	
	if (is.null(annotationLibrary)) {
		annotationLibrary <- annotation(eset)}
	
	annotationPkg <- paste(annotation(eset), ".db", sep="")
	require(annotationPkg, character.only = TRUE)
	fNames <- featureNames(eset)
	
	### ENTREZID
	pData(featureData(eset))[,1] <- unlist(AnnotationDbi:::mget(fNames, eval(parse(text=paste(annotationLibrary, "ENTREZID", sep="")))))[fNames]   
	colnames(pData(featureData(eset)))[1] <- "ENTREZID"
	fvarMetadata(eset)[1,1] <- "Entrez ID as retrieved from annotation package"
	### ENSEMBL ID
	pData(featureData(eset))[,2] <- unlist(AnnotationDbi:::mget(fNames,eval(parse(text=paste(annotationLibrary, "ENSEMBL", sep="")))))[fNames]
	colnames(pData(featureData(eset)))[2] <- "ENSEMBLID"
	fvarMetadata(eset)[2,1] <- "Ensembl ID as retrieved from annotation package"
	### Gene Symbol
	pData(featureData(eset))[,3] <- unlist(AnnotationDbi:::mget(fNames,eval(parse(text=paste(annotationLibrary, "SYMBOL", sep="")))))[fNames]
	colnames(pData(featureData(eset)))[3] <- "SYMBOL"
	fvarMetadata(eset)[3,1] <- "Gene symbol as retrieved from annotation package"
	### Description
	pData(featureData(eset))[,4] <- unlist(AnnotationDbi:::mget(fNames,eval(parse(text=paste(annotationLibrary, "GENENAME", sep="")))))[fNames]
	colnames(pData(featureData(eset)))[4] <- "GENENAME"
	fvarMetadata(eset)[4,1] <- "Description as retrieved from annotation package"
	return(eset)
}

