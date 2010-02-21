\name{spectralMap}
\alias{spectralMap}
\title{Draw a Spectral Map}
\description{
  Generic function to draw a spectral map.
}
\usage{
spectralMap(object, groups, ...)
}
\arguments{
  \item{object}{object of class ExpressionSet}
  \item{groups}{string indicating the name of the column in the phenoData that
     defines the groups}
  \item{\dots}{further arguments to be passed to the methods}
}
\value{
  Object of class \code{plot.mpm}, i.e. the S3 output object of the \code{plot.mpm}
  function of the \code{mpm} package
}
\note{Coloring of groups on the spectralMap uses the a4 palette as produced
      by \code{a4palette}}
\references{Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U., Molenberghs, G.,
  Lewi, P.J. (2003). Graphical exploration of gene expression data: a
  comparative study of three multivariate methods. \emph{Biometrics}
  \bold{59}, 1131-1140.
  
  Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
    Microarrays, Chapman \& Hall/CRC, pp. 148 - 153.
}
\author{Tobias Verbeke}
\seealso{\code{\link{spectralMap-methods}}, \code{\link[mpm]{plot.mpm}}}

\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)

  spectralMap(object = ALL, groups = "BT", legendPos = 'bottomright')

  spectralMap(object = ALL, groups = "BT",
      plot.mpm.args = list(label.tol = 10, rot = c(-1, 1), sub = "", lab.size = 0.65,
	    dim = c(1,2), sampleNames = FALSE, zoom = c(1,5), col.size = 2, 
	    do.smoothScatter = TRUE))
	  
  spectralMap(object = ALL, groups = "BT",
	  plot.mpm.args = list(label.tol = 10, rot = c(-1, 1), sub = "", lab.size = 0.65,
	      dim = c(1,2), sampleNames = as.character(pData(ALL)$BT),
          zoom = c(1,5), col.size = 2, do.smoothScatter = TRUE))
}
}
\keyword{hplot}
%\item{makeLognormal}{boolean indicating whether one wants to exponentiate the
%     data to make them lognormally shaped (\code{TRUE}; the default) or not 
%     (\code{FALSE})}
%  \item{plot.mpm.args}{list of arguments that can be passed to the 
%     \code{plot.mpm} function that actually draws the plot}
%  \item{probe2gene}{boolean indicating whether one wants to display the gene symbols
%     for the labeled points (\code{TRUE}) or not (\code{FALSE}; the default) }
  