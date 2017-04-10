\name{SDOwave}
\alias{SDOwave}
\docType{data}
\title{
Longitudinal Social Dominance Orientation (SDO)
}
\description{
Contains 4 SDO items measured across 5 years (1996-2000). Each item is scored on a 7-point scale. 
}
\usage{
data("SDOwave")
}
\format{
Data frame containing 612 subjects, 4 items measure across 5 years (wide format). Here are the item labels for one year:
\describe{
	\item{\code{I1.1996}}{It's probably a good thing that certain groups are at the top and other groups are
at the bottom.}
  \item{\code{I2.1996}}{Inferior groups should stay in their place.}
  \item{\code{I3.1996}}{We should do what we can to equalize conditions for different groups (reversed).}
  \item{\code{I4.1996}}{Increased social equality is beneficial to society (reversed).}
}
}

\references{
Sidanius, J., & Pratto, F. (2001). Social Dominance: An Intergroup Theory of Social Hierarchy and Oppression. Cambridge University Press, Cambridge, UK.

}

\note{
Thanks to Jim Sidanius for providing this dataset.
}



\examples{
data("SDOwave")
str(SDOwave)
}
\keyword{datasets}
