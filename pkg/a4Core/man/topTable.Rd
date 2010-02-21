\name{topTable}
\alias{topTable}
\title{
  S4 Generic for obtaining a top table 
}
\description{
  a top table is a rectangular object (e.g. data frame) which
  lists the top n most relevant variables
}
\usage{
topTable(fit, n, ...)
}
\arguments{
  \item{fit}{
    object for which to obtain a top table, generally a fit object for a given model class
}
  \item{n}{
    number of features (variables) to list in the top table, ranked by importance
}
  \item{\dots}{
    further arguments for specific methods
}
}
\author{
  Tobias Verbeke
}
\keyword{models}
