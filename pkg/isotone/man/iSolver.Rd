\name{iSolver}
\alias{iSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{SILF Loss}
\description{Minimizes soft insensitive loss function (SILF) for support vector regression.
}
\usage{
iSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector, \code{weights} 
  with optional observation weights, \code{beta} between 0 and 1, and \code{eps} > 0}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = iSolver}.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\references{Efron, B. (1991). Regression percentiles using asymmetric squared error loss. Statistica Sinica, 1, 93-125.}

\seealso{\code{\link{activeSet}}}
\examples{

##Fitting isotone regression using active set
set.seed(12345)
y <- rnorm(9)              ##response values
w <- rep(1,9)              ##unit weights
eps <- 2
beta <- 0.4

btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.silf <- activeSet(btota, iSolver, weights = w, y = y, beta = beta, eps = eps)

}
\keyword{models}
