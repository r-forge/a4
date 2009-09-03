\name{ExpressionSetWithComputation-class}
\Rdversion{1.1}
\docType{class}
\alias{ExpressionSetWithComputation-class}

\title{Class "ExpressionSetWithComputation"}
\description{This class adds statistical information to the exprs of the ExpressionSet
  as well as descriptive information to the pData of the ExpressionSet}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ExpressionSetWithComputation", assayData, phenoData, featureData, experimentData, annotation, exprs, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}}
    \item{\code{annotation}:}{Object of class \code{"character"}}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{\link[Biobase]{ExpressionSet}}, directly.
Class \code{\link[Biobase]{eSet}}, by class "ExpressionSet", distance 2.
Class \code{\link[Biobase]{VersionedBiobase}}, by class "ExpressionSet", distance 3.
Class \code{\link[Biobase]{Versioned}}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
No methods defined with class "ExpressionSetWithComputation" in the signature.
}
\author{Tobias Verbeke}
\seealso{
	\code{\link[Biobase]{ExpressionSet}}, \code{\link{computeDEStats}}
}
\keyword{classes}
