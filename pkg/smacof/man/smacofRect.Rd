\name{smacofRect}
\alias{smacofRect}

\title{Rectangular smacof}
\description{Variant of smacof for rectangular matrices (typically ratings, preferences) which is also known as metric unfolding.
}
\usage{
smacofRect(diss, ndim = 2, weightmat = NULL, init = NULL, verbose = FALSE,
           itmax = 100, reg = 1e-6, eps = 1e-6)
}

\arguments{
  \item{diss}{Data frame or matrix of preferences, ratings, dissimilarities.}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{itmax}{Maximum number of iterations}
  \item{reg}{Scalar for configuration distances}
  \item{eps}{Convergence criterion}
}
\details{Creates an object of class \code{smacofR}. 

}

\value{
  \item{confdiss}{Configuration dissimilarities}
  \item{conf.row}{Matrix of final row configurations}
  \item{conf.col}{Matrix of final column configurations}
  \item{stress}{Final stress value}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nind}{Number of individuals (rows)}
  \item{nobj}{Number of objects (columns)}
}
\references{de Leeuw, J. \& Mair, P. (2008). Multidimensional scaling using majorization: The R package smacof.}
\author{Jan de Leeuw and Patrick Mair}

\seealso{to be added}

\examples{
data(breakfast)
res <- smacofRect(breakfast)
res
summary(res)
}

\keyword{models}