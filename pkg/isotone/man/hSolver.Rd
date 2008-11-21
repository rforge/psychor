\name{hSolver}
\alias{hSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Huber Loss Function}
\description{Solver for Huber's robust loss function.
}
\usage{
hSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector and \code{weights} 
  with optional observation weights, and \code{eps}}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = hSolver}.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\references{
Huber, P. (1982). Robust Statistics. Chichester: Wiley.
}

\seealso{\code{\link{activeSet}}}
\examples{

##Fitting isotone regression using active set
set.seed(12345)
z <- 9:1                   ##predictor values
y <- rnorm(9)              ##response values
w <- rep(1,9)              ##unit weights
eps = 0.01
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.huber <- activeSet(z, btota, hSolver, weights = w, y = y, eps = eps)

}
\keyword{models}
