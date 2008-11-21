\name{oSolver}
\alias{oSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lp norm}
\description{Solver for Lp-norm.
}
\usage{
oSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector, \code{weights} 
  as an optional weight vector, and \code{p} as the exponent for the Lp-norm.}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = oSolver}.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\seealso{\code{\link{activeSet}}}
\examples{

##Fitting isotone regression 
set.seed(12345)
z <- 9:1                   ##predictor values
y <- rnorm(9)              ##response values
w <- diag(rep(1,9))        ##unit weight matrix
p <- 3
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
#fit.pnorm <- activeSet(z, btota, oSolver, weights = w, y = y, p = p)

}
\keyword{models}
