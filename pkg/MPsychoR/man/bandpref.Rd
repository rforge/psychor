\name{bandpref}
\alias{bandpref}
\docType{data}
\title{
Band Preferences
}
\description{
Toy dataset involving paired comparisons of bands. 200 people stated their preferences of 5 bands in a paired comparison design (no undecided answer allowed). 
}

\usage{
data("bandpref")
}

\format{
A data frame with 10 paired comparisons (200 people): 
  \describe{
	\item{\code{Band1}}{First band}
	\item{\code{Band2}}{Second band}
	\item{\code{Win1}}{How often first band was preferred}
	\item{\code{Win2}}{How often second band was preferred}
	}
}


\examples{
data("bandpref")
str(bandpref)
}
\keyword{datasets}
