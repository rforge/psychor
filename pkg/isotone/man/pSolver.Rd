\name{pSolver}
\alias{pSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantile Regression}
\description{Solver for the general p-quantile monotone regression problem with optional weights.
}
\usage{
pSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector, \code{weights} 
  with optional observation weights, \code{aw} and \code{bw} as quantile weights.}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = pSolver}. Note
that if \code{aw} = \code{bw}, we get the weighted median and therefore we solved the weighted absolute norm.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\references{Koenker, R. (2005). Quantile regression. Cambridge, MA: Cambridge University Press.
}

\seealso{\code{\link{activeSet}}}
\examples{

##Fitting quantile regression
set.seed(12345)
y <- rnorm(9)              ##response values
w <- rep(1,9)              ##unit weights
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.p <- activeSet(btota, pSolver, weights = w, y = y, aw = 0.3, bw = 0.7)

}
\keyword{models}
