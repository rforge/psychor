\name{gpava}
\alias{gpava}
\alias{print.pava}
\alias{weighted.median}
\alias{weighted.fractile}
\alias{plot.pava}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Pooled-Adjacent-Violators Algorithm (PAVA)}
\description{Pooled-adjacent-violators algorithm for general isotone regression problems. It allows for general convex target function, 
multiple measurements, and different approaches for handling ties.}
\usage{
gpava(z, y, weights = NULL, solver = weighted.mean, ties = "primary", a = NA, b = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector of abscissae values}
  \item{y}{Vector or list of vectors of responses}
  \item{weights}{Vector of list of vectors of observation weights}
  \item{solver}{Either \code{weighted.mean}, \code{weighted.median}, \code{weighted.fractile}, or 
  a user-specified function (see below)}
  \item{ties}{Treatment of ties, either "primary", "secondary", or "tertiary"}
  \item{a}{Fractile specification for \code{weighted.fractile} solver}
  \item{b}{Fractile specification for \code{weighted.fractile} solver}
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
  \item{x}{Fitted values}
  \item{y}{Observed response}
  \item{z}{Observed predictors}
  \item{w}{Weights}
  \item{solver}{Convex function}
  \item{call}{Matched call}
}
\references{
de Leeuw, J., Hornik, K., Mair, P. (2008). Isotone Optimization in R: Pool-Adjacent-Violators 
Algorithm (PAVA) and Active Set Methods. Journal of Statistical Software, forthcoming. 
}
\author{Kurt Hornik, Jan de Leeuw, Patrick Mair}


\examples{

data(pituitary)
##different tie approaches
gpava(pituitary[,1],pituitary[,2], ties = "primary")
gpava(pituitary[,1],pituitary[,2], ties = "secondary")
gpava(pituitary[,1],pituitary[,2], ties = "tertiary")

##different target functions
gpava(pituitary[,1],pituitary[,2], solver = weighted.mean)
gpava(pituitary[,1],pituitary[,2], solver = weighted.median)
gpava(pituitary[,1],pituitary[,2], solver = weighted.fractile, a = 0.25, b = 0.75)


##repeated measures
data(posturo)
res <- gpava(posturo[,1],posturo[,2:4], ties = "secondary")
plot(res)

}
\keyword{models}
