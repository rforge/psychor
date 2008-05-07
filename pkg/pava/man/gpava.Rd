\name{gpava}
\alias{gpava}
\alias{print.pava}
\alias{weighted.median}
\alias{weighted.pth.fractile}
\alias{plot.pava}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized PAVA}
\description{PAVA for general convex functions, including 
}
\usage{
gpava(x, y, w = NULL, solver = weighted.mean, ties = "primary")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of list of vectors of predictors}
  \item{y}{Vector of list of vectors of responses}
  \item{w}{Vector of list of vectors of observation weights}
  \item{solver}{Either \code{weighted.mean}, \code{weighted.median}, \code{weighted.pth.fractile}, or 
  a user-specified function (see below)}
  \item{ties}{Treatment of ties, either "primary", "secondary", or "tertiary"}
}
\details{A Pool Adjacent Violators Algorithm framework for minimizing problems
like 

\deqn{\sum_i \sum_{J_i} w_{ij} f(y_{ij}, m_i)}

under the constraint \eqn{m_1 \le ... \le m_n} with \eqn{f} a convex function in m. Note that this formulation allows for repeated data in each block
(i.e. each list element of \code{y}, and hence is more general than the usual pava/isoreg ones.  

A solver for the unconstrained \eqn{\sum_k w_k f(y_k, m) -> min!} can be specified. 
Typical cases are \eqn{f(y, m) = |y - m|^p} for \eqn{p = 2} (solved by weighted mean) and \eqn{p = 1} (solved by weighted median), respectively. 
The user can also specify his own function \code{foo(y, w)} with responses and weights as arguments. It 
should return a single numerical value. 
}

\value{
  Generates an object of class \code{gpava}.
  \item{yfit}{Fitted values}
  \item{y}{Observed response}
  \item{x}{Observed predictors}
  \item{w}{Weights}
  \item{solver}{Convex function}
  \item{call}{Matched call}
}
\references{
de Leeuw, J., Hornik, K., Mair, P. (2008). Pool-Adjacent-Violators Algorithm (PAVA) in R. 
Forthcoming.
}
\author{Kurt Hornik, Jan de Leeuw, Patrick Mair}


\examples{

data(pituitary)
res1 <- gpava(pituitary[,1],pituitary[,2], ties = "primary")
res2 <- gpava(pituitary[,1],pituitary[,2], ties = "secondary")
res3 <- gpava(pituitary[,1],pituitary[,2], ties = "tertiary")

##now with three different functions
res.mean <- gpava(pituitary[,1],pituitary[,2], solver = weighted.mean)
res.median <- gpava(pituitary[,1],pituitary[,2], solver = weighted.median)

##repeated measures
data(posturo)
res <- gpava(posturo[,1],posturo[,2:4], ties = "secondary")
plot(res)

}
\keyword{models}
