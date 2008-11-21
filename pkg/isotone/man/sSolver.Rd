\name{sSolver}
\alias{sSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Negative Poisson Log-Likelihood}
\description{Solver for the negative Poisson log-likelihood
}
\usage{
sSolver(z, a, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed response}
  \item{a}{Matrix with active constraints}
  \item{extra}{List with element \code{y} containing the observed response vector}
}
\details{This function is called internally in \code{activeSet} by setting \code{mySolver = sSolver}.
}

\value{
  \item{x}{Vector containing the fitted values}
  \item{lbd}{Vector with Lagrange multipliers}
  \item{f}{Value of the target function}
  \item{gx}{Gradient at point x}
}

\seealso{\code{\link{activeSet}}}
\examples{

##Minimizing Poisson log-liklihood
set.seed(123)
z <- 9:1                   ##predictor values
y <- sample(1:100, 9)      ##response values
btota <- cbind(1:8, 2:9)   ##Matrix defining isotonicity (total order)
fit.poiss <- activeSet(z, btota, sSolver, y = y)
}
\keyword{models}
