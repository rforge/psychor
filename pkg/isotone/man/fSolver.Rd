\name{fSolver}
\alias{fSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{User-Specified Loss Function}
\description{Specification of a differentiable convex loss function. 
}
\usage{
fSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{fobj} containing the target function and \code{gobj} with the 
  first derivative}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = fSolver}. It uses
\code{optim()} with \code{"BFGS"} for optimization.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\seealso{\code{\link{activeSet}}}
\examples{

##Fitting isotone regression using active set (L2-norm user-specified)
set.seed(12345)
y <- rnorm(9)              ##response values
w <- rep(1,9)              ##unit weights
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.convex <- activeSet(btota, fSolver, fobj = function(x) sum(w*(x-y)^2), 
gobj = function(x) 2*drop(w*(x-y)), y = y, weights = w)
}
\keyword{models}
