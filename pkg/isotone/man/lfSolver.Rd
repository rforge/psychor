\name{lfSolver}
\alias{lfSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{General Least Squares Loss Function}
\description{Solver for the general least squares monotone regression problem of the form
(y-x)'W(y-x). 
}
\usage{
lfSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector and \code{weights} 
  as weight matrix W which is not necessarily positive definite.}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = lfSolver}.
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
y <- rnorm(9)              ##response values
w <- diag(rep(1,9))        ##unit weight matrix
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
#fit.lf <- activeSet(btota, lfSolver, weights = w, y = y)

}
\keyword{models}
