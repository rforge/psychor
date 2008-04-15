\name{gpava}
\alias{gpava}
\alias{print.pava}
\alias{weighted.median}
\alias{weighted.pth.fractile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized PAVA}
\description{PAVA for general convex functions, including 
}
\usage{
gpava(x, y, w = NULL, solver = weighted.mean, merger = c, ties = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of list of vectors of predictors}
  \item{y}{Vector of list of vectors of responses}
  \item{w}{Vector of list of vectors of observation weights}
  \item{solver}{Either \code{weighted.mean}, \code{weighted.median}, \code{weighted.pth.fractile}, or 
  a user-specified function (see below)}
  \item{merger}{Function for merging data}
  \item{ties}{Treatment of ties, either "none", "primary", "secondary", or "tertiary"}
}
\details{A Pool Adjacent Violators Algorithm framework for minimizing problems
like 

\deqn{\sum_i \sum_{J_i} w_{ij} f(y_{ij}, m_i)}

under the constraint \eqn{m_1 < ... < m_n} with \eqn{f} a convex function in m. Note that this formulation allows for repeated data in each block
(i.e. each list element of \code{y}, and hence is more general than the usual pava/isoreg ones.  

A solver for the unconstrained \eqn{\sum_k w_k f(y_k, m) -> min!} is needed. Typical cases are \eqn{f(y, m) = |y - m|^p} for \eqn{p = 2} 
(solved by weighted mean) and \eqn{p = 1} (solved by weighted median), respectively. The user 
can specify his own function \code{foo(y, w)} with responses and weights as arguments. It 
should return a single numerical value. 
}

\value{
  Generates an object of class \code{pava}.
  \item{yfit}{Fitted values}
  \item{y}{Observed response}
  \item{x}{Observed predictors}
}
\references{
de Leeuw, J., Hornik, K., Mair, P. (2008). Pool-Adjacent-Violators Algorithm (PAVA) in R. 
Forthcoming.
}
\author{Kurt Hornik, Jan de Leeuw, Patrick Mair}
\examples{

x <- c(1,1,2,2,2,3,3,3,3)
y <- c(1,4,0,3,3,5,6,7,4)

gpava(x, y)                          ##no tie (same result as isoreg()
gpava(x, y, ties = "primary")        ##primary approach
gpava(x, y, ties = "secondary")      ##secondary approach
gpava(x, y, ties = "tertiary")       ##tertiary appraoch

##now with three different functions
gpava(1:9, y, ties = "primary", solver = weighted.mean)
gpava(1:9, y, ties = "primary", solver = weighted.median)
gpava(1:9, y, ties = "primary", solver = weighted.pth.fractile)

}
\keyword{models}
