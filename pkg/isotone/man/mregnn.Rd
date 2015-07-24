\name{mregnn}
\alias{mregnn}
\alias{mregnnM}
\alias{mregnnP}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression with Linear Inequality Restrictions on Predicted Values}
\description{The package contains three functions for fitting regressions with inequality restrictions:
\code{mregnn} is the most general one, allowing basically for any partial orders, \code{mregnnM} poses a monotone restriction on the fitted values, \code{mregnnP} restricts the predicted values to be positive. Monre details can be found below. 
}
\usage{
mregnn(x, y, a)
mregnnM(x, y)
mregnnP(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Can be a spline basis.}
  \item{y}{Response.}
  \item{a}{Matrix containing order restrictions.}
}
\details{These functions solve the problem 
\deqn{
f(b) = \frac{1}{2}(y - Xb)'(y - Xb)
}
over all \eqn{b} for which \eqn{A'Xb \geq 0}. \eqn{A} can be used require the transformation to be non-negative, or increasing, or satisfying any partial order. 
}

\value{
  \item{xb}{Predicted values.}
  \item{lb}{Solution of the dual problem.}
  \item{f}{Value of the target function}
}
\references{
de Leeuw, J. (2015). Regression with Linear Inequality Restrictions on Predicted Values. \url{http://rpubs.com/deleeuw/78897}.
}
\examples{
## Compute the best fitting quadratic polynomial (in black) 
## and monotone quadratic polynomial (in blue)
set.seed(12345)
x <- outer(1:10,1:3,"^")
x <- apply(x,2,function(x)
x - mean(x))
x <- apply (x,2,function(x)
x / sqrt (sum(x ^ 2)))
y <- rowSums(x) + rnorm(10)
plot(x[,1], y, lwd = 3, col = "RED", xlab = "x", ylab = "P(x)")
o <- mregnnM(x,y)
lines(x[,1], o$xb, col = "BLUE", lwd = 2)
xb <- drop(x \%*\% qr.solve(x,y))
lines(x[,1],xb,col="BLACK", lwd = 2)

## same monotone model through basic mregnn()
difmat <- function (n) {
  m1 <- ifelse(outer(1:(n - 1),1:n,"-") == -1, 1, 0)
  m2 <- ifelse(outer(1:(n - 1),1:n,"-") == 0,-1, 0)
  return (m1 + m2)
}
a <- difmat(nrow(x))      ## order restriction
o2 <- mregnn(x, y, a)
}
\keyword{models}
