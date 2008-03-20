\name{smacofConstraint}
\alias{smacofConstraint}

\title{SMACOF Constraint}
\description{
SMACOF with constraints on the configuration
}
\usage{
smacofConstraint(delta, constraint = "linear", external, ndim = 2, weightmat = NULL, 
metric = TRUE, ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{constraint}{Type of constraint: \code{"linear"}, \code{"unique"}, \code{"diagonal"}, or a user-specified function (see details)}
  \item{external}{Data frame or matrix with external covariates}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{metric}{If \code{FALSE} non-metric MDS is performed}
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}
\details{The user can specify a function with the following arguments: configuration matrix
(n x p), matrix V (based on the weight matrix, see vignette), external scale matrix. 
The function must return a (n x p) matrix of resulting configurations. In the examples we show how linearly 
restricted configurations can be specified by the user.  

}
\value{
  \item{obsdiss}{Observed dissimilarities, normalized}
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

\seealso{\code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}
\examples{

## SMACOF with linear configuration constraints
data(kinshipdelta)
data(kinshipscales)
res.lin1 <- smacofConstraint(kinshipdelta, constraint = "linear", external = kinshipscales)

## The same model, linear restrictions specified externally
## constrfun <- function(X, V, external) {
##      return(external%*%solve(crossprod(external,V%*%external),crossprod(external,V%*%X)))
## }
## res.lin2 <- smacofConstraint(kinshipdelta, constraint = constrfun, external = kinshipscales)

}

\keyword{models}
