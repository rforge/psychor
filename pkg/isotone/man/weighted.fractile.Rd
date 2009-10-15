\name{weighted.fractile}
\alias{weighted.fractile}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Median}
\description{Computes the weighted fractile of a numeric vector
}
\usage{
weighted.fractile(y, w, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A numeric vector containing the values whose fractile is to be computed}
  \item{w}{A vector of length \code{y} giving the weights to use for each element of \code{y}}
  \item{p}{Fractile specification; value between 0 and 1}
}

\seealso{\code{\link{weighted.mean}}, \code{\link{weighted.median}}}
\examples{
y <- 1:9
w <- c(rep(1,5), rep(2,4))
res <- weighted.fractile(y, w, p = 0.33)
}
\keyword{models}
