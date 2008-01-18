\name{smacofSym}
\alias{smacofSym}

\title{Symmetric smacof}
\description{
Basic smacof on symmetric dissimilarity matrix 
}
\usage{
smacofSym(diss, ndim = 2, weightmat = NULL, init = NULL, metric = TRUE, ties = "primary", verbose = FALSE, relax = 1, modulus = 1, itmax = 100, eps = 1e-06)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diss}{Either a symmetric dissimilarity matrix or an object of class \code{dist}}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{metric}{If \code{TRUE} metric MDS is performed, if \code{FALSE} non-metric MDS}
  \item{ties}{Tie specification for PAVA (for non-metric MDS only)}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{Block relaxation}
  \item{modulus}{Modulus for non-metric MDS}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}
\details{
to be added
}
\value{
  \item{obsdiss}{Observed dissimilarities}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{stress.m}{stress value for metric MDS}
  \item{stress.nm}{stress value for non-metric MDS (if computed)}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. \& Mair, P. (2008). Multidimensional scaling using majorization: The R package smacof.}
\author{Jan de Leeuw and Patrick Mair}

\seealso{to be added}
\examples{
data(ekman)
res <- smacofSym(ekman)
res
summary(res)
}

\keyword{models}
