\name{topTable-methods}
\docType{methods}
\alias{topTable-methods}
\alias{topTable,glmnet-method}
\title{Methods for topTable}
\description{
  Methods for topTable. topTable extracts the top n most important features
  for a given classification or regression procedure 
}
\section{Methods}{
\describe{

glmnet

\item{fit = "glmnet", n = "numeric"}{glmnet objects are produced by \code{lassoClass} or \code{lassoReg}}

}
}
\arguments{
  \item{fit}{object resulting from a classification or regression procedure}
  \item{n}{number of features that one wants to extract from a table that
    ranks all features according to their importance in the classification
    or regression model; defaults to 10 for limma objects}
}
\keyword{methods}
\keyword{manip}

