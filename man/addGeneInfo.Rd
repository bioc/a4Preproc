% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addGeneInfo.R
\name{addGeneInfo}
\alias{addGeneInfo}
\title{Utility Function to Add Annotation to existing ExpressionSet Objects}
\usage{
addGeneInfo(eset, annotationLibrary = NULL)
}
\arguments{
\item{eset}{ExpressionSet object for to which one wants to add 
additional annotation information}

\item{annotationLibrary}{Annotation Library to use. Must be specified when working with Entrez gene IDs.
In this case, one can make use of the JnJ annotation packages such as \code{hgu133plus2hsentrezgJnJ}.
If not specified, the annotation of the package will be automatically requested
with \code{annotation()} of the expressionSet object \code{eset} and then Affymetrix probe set IDs are expected in \code{featureNames}}
}
\value{
a new ExpressionSet object with the additional information stored as feature data
}
\description{
Adds appropriate featureData, similar to the metadata added in the pipeline script,
to the ExpressionSet object.
}
\details{
Slots of featureData(a4ALL) are
\itemize{
\item \kbd{Entrez ID}~: Entrez ID as retrieved from annotation package
\item \kbd{Ensembl ID}~: Ensembl ID as retrieved from annotation package
\item \kbd{Gene Symbol}~: Gene symbol as retrieved from annotation package
\item \kbd{Description}~: Description as retrieved from annotation package 
}
}
\note{
One should always use subscripting of featureData by column name
(e.g. \code{featureData(a4ALL)$`Entrez ID`}; as the pipeline 
ExpressionSets have one additional column compared to the ExpressionSet
objects produced by \code{addGeneInfo}, i.e. column 2 of the pipeline ExpressionSets
corresponds to column one of an \code{addGeneInfo} ExpressionSet.
}
\examples{
  
library(ALL)
data(ALL)
a4ALL <- addGeneInfo(ALL)
head(featureData(a4ALL)$`Entrez ID`)
}
\author{
Tobias Verbeke, Steven Osselaer
}
\keyword{manip}
