\name{aSolver}
\alias{aSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Asymmetric Least Squares}
\description{Minimizes Efron's asymmetric least squares regression.
}
\usage{
aSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector, \code{weights} 
  with optional observation weights, weight \code{aw} for y > x, and weight \code{bw} for y <= x}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = aSolver}.
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
z <- 9:1                   ##predictor values
y <- rnorm(9)              ##response values
w <- rep(1,9)              ##unit weights
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.asy <- activeSet(z, btota, aSolver, weights = w, y = y, aw = 0.3, bw = 0.5)

}
\keyword{models}
