\name{knotsQ}
\alias{knotsQ}
\alias{knotsR}
\alias{knotsD}
\alias{knotsE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility functions for knot specifications. 
}
\description{
There is \code{knotsQ} for knots at the quantiles, \code{knotsR} for knots equally spaced on the range, \code{knotsD} for knots at the data points, and \code{knotsE} for no interior knots.
}
\usage{
knotsQ(x, n = 5) 
knotsR(x, n = rep(5, ncol (x)))
knotsE(x)
knotsD(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data frame}
  \item{n}{Number of knots}
}
\details{
Simple linear or polynomial transformations can be specified using \code{knotsE} and the corresponding degree in the subsequent Gifi function call. If the data frame contains factors, they are converted internally to numeric 
using \code{makeNumeric}. 
}
\value{
A list containing knots for each variable.
}

\examples{
ABC6 <- ABC[,6:11]
knotsQ(ABC6)        ## tercile knots
knotsQ(ABC6, n = 6) ## quartile knots 
knotsR(ABC6)        ## equally spaced knots
knotsE(ABC6)        ## no interior knots
knotsD(ABC6)        ## knots at the data points
}
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")
