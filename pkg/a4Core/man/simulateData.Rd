\name{simulateData}
\alias{simulateData}
\title{Simulate Data for Package Testing and Demonstration Purposes}
\description{
  Simulate Data for Package Testing and Demonstration Purposes}
\usage{
  simulateData(nCols = 40, nRows = 1000, nEffectRows = 5, nNoEffectCols = 5,
      betweenClassDifference = 1, withinClassSd = 0.5)
}
\arguments{
  \item{nCols}{number of samples; currently this should be an even number}
  \item{nRows}{number of features (genes)}
  \item{nEffectRows}{number of differentially expressed features}
  \item{nNoEffectCols}{number of samples for which the profile
    of a differentially expressed feature will be set similar 
    to the other class} 
  \item{betweenClassDifference}{Average mean difference between the two classes 
    to simulate a certain signal in the features for which an effect was introduced;
    the default is set to 1}
  \item{withinClassSd}{Within class standard deviation used to add a certain noise 
    level to the features for which an effect was introduced; the default standard
    deviation is set to 0.5}
}
\value{
  object of class ExpressionSet with the characteristics specified 
}
\note{
  The simulation assumes the variances are equal between the two classes.
  Heterogeneity could easily be introduced in the simulation if this would
  be requested by the users.
}
\author{W. Talloen and T. Verbeke}
\examples{
  someEset <- simulateData(nCols = 40, nRows = 1000, nEffectRows = 5, nNoEffectCols = 5)
  someEset
}
\keyword{manip}

