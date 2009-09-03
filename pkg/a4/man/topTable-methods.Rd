\name{topTable-methods}
\docType{methods}
\alias{topTable-methods}
\alias{topTable,glmnet-method}
\alias{topTable,limma-method}
\alias{topTable,MArrayLM-method}
\alias{topTable,pamClass-method}
\alias{topTable,rfClass-method}
\alias{topTable,tTest-method}
\title{Methods for topTable}
\description{
  Methods for topTable. topTable extracts the top n most important features
  for a given classification or regression procedure 
}
\section{Methods}{
\describe{

glmnet

\item{fit = "glmnet", n = "numeric"}{glmnet objects are produced by \code{lassoClass} or \code{lassoReg}}

limma

\item{fit = "limma", n = "numeric"}{limma objects are produced by \code{limma2Groups}}

MarrayLM

\item{fit = "limma", n = "numeric"}{MarrayLM objects are produced by \code{lmFit} of the \code{limma package}}


pamClass

\item{fit = "pamClass", n = "numeric"}{pamClass objects are produced by \code{pamClass}}

rfClass

\item{fit = "rfClass", n = "numeric"}{rfClass objects are produced by \code{rfClass}}

tTest

\item{fit = "tTest", n = "numeric"}{tTest objects are produced by \code{tTest}}


}
}
\arguments{
  \item{fit}{object resulting from a classification or regression procedure}
  \item{n}{number of features that one wants to extract from a table that
    ranks all features according to their importance in the classification
    or regression model}
}
\keyword{methods}
\keyword{manip}
