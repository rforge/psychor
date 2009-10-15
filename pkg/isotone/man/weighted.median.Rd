\name{weighted.median}
\alias{weighted.median}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Median}
\description{Computes a weighted median of a numeric vector
}
\usage{
weighted.median(y, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A numeric vector containing the values whose median is to be computed}
  \item{w}{A vector of length \code{y} giving the weights to use for each element of \code{y}}
}

\seealso{\code{\link{weighted.mean}}, \code{\link{weighted.fractile}}}
\examples{
y <- 1:9
w <- c(rep(1,5), rep(2,4))
res <- weighted.median(y, w)
}
\keyword{models}
