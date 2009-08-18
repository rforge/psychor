\name{dSolver}
\alias{dSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Absolute Value Norm}
\description{Solver for the least absolute value norm with optional weights.
}
\usage{
dSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector and \code{weights} 
  with optional observation weights}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = dSolver}.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\seealso{\code{\link{activeSet}}}
\examples{

##Fitting weighted absolute norm problem
set.seed(12345)
y <- rnorm(9)              ##response values
w <- rep(1,9)              ##unit weights
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.abs <- activeSet(btota, dSolver, weights = w, y = y)

}
\keyword{models}
