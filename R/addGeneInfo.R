

### utility function to transform an ExpressionSet into an ExpressionSet with meta data similar
### to the ExpressionSets used in the pipeline

#' Utility Function to Add Annotation to existing ExpressionSet Objects
#' 
#' Adds appropriate featureData, similar to the metadata added in the pipeline script,
#' to the ExpressionSet object.
#' @param eset ExpressionSet object for to which one wants to add 
#' additional annotation information
#' @param annotationLibrary Annotation Library to use. Must be specified when working with Entrez gene IDs.
#'   In this case, one can make use of the JnJ annotation packages such as \code{hgu133plus2hsentrezgJnJ}.
#'   If not specified, the annotation of the package will be automatically requested
#'   with \code{annotation()} of the expressionSet object \code{eset} and then Affymetrix probe set IDs are expected in \code{featureNames}
#' @details
#' Slots of featureData(a4ALL) are
#' \itemize{
#' \item \kbd{Entrez ID}~: Entrez ID as retrieved from annotation package
#' \item \kbd{Ensembl ID}~: Ensembl ID as retrieved from annotation package
#' \item \kbd{Gene Symbol}~: Gene symbol as retrieved from annotation package
#' \item \kbd{Description}~: Description as retrieved from annotation package 
#' }
#' @return a new ExpressionSet object with the additional information stored as feature data
#' @note One should always use subscripting of featureData by column name
#' (e.g. \code{featureData(a4ALL)$`Entrez ID`}; as the pipeline 
#' ExpressionSets have one additional column compared to the ExpressionSet
#' objects produced by \code{addGeneInfo}, i.e. column 2 of the pipeline ExpressionSets
#' corresponds to column one of an \code{addGeneInfo} ExpressionSet.
#' @author Tobias Verbeke, Steven Osselaer
#' @examples   
#' library(ALL)
#' data(ALL)
#' a4ALL <- addGeneInfo(ALL)
#' head(featureData(a4ALL)$`Entrez ID`)
#' @keywords manip
#' @importFrom BiocGenerics annotation mget
#' @importFrom Biobase featureNames featureData `featureData<-` `pData<-` pData `fvarMetadata<-` fvarMetadata
#' @export
addGeneInfo <- function(eset, annotationLibrary = NULL){
	
	if (is.null(annotationLibrary)) {
		annotationLibrary <- annotation(eset)
	}
	
	annotationPkg <- paste(annotationLibrary, ".db", sep="")
	require(annotationPkg, character.only = TRUE)
	fNames <- featureNames(eset)
	
	### ENTREZID
	pData(featureData(eset))[,1] <- unlist(mget(fNames, eval(parse(text=paste(annotationLibrary, "ENTREZID", sep="")))))[fNames]   
	colnames(pData(featureData(eset)))[1] <- "ENTREZID"
	fvarMetadata(eset)[1,1] <- "Entrez ID as retrieved from annotation package"
	### ENSEMBL ID
	pData(featureData(eset))[,2] <- unlist(mget(fNames,eval(parse(text=paste(annotationLibrary, "ENSEMBL", sep="")))))[fNames]
	colnames(pData(featureData(eset)))[2] <- "ENSEMBLID"
	fvarMetadata(eset)[2,1] <- "Ensembl ID as retrieved from annotation package"
	### Gene Symbol
	pData(featureData(eset))[,3] <- unlist(mget(fNames,eval(parse(text=paste(annotationLibrary, "SYMBOL", sep="")))))[fNames]
	colnames(pData(featureData(eset)))[3] <- "SYMBOL"
	fvarMetadata(eset)[3,1] <- "Gene symbol as retrieved from annotation package"
	### Description
	pData(featureData(eset))[,4] <- unlist(mget(fNames,eval(parse(text=paste(annotationLibrary, "GENENAME", sep="")))))[fNames]
	colnames(pData(featureData(eset)))[4] <- "GENENAME"
	fvarMetadata(eset)[4,1] <- "Description as retrieved from annotation package"
	return(eset)
}

