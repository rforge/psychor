\name{monspline}
\alias{monspline}
\alias{bsb}
\title{Monotone Splines}
\description{
  Here we provide the following two utiliy function: \code{bsb()} can be used to create a B-spline basis. \code{monspline()} transforms an arbitrary spline (non-monotone) into a monotone spline by using a dual method that uses \code{nnls()} to find the Lagrange multipliers 
  (see description for details). 
}
\usage{
bsb(x, knots, order = 3) 

monspline(y, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector containing the predictor variable.}
  \item{knots}{A vector specifying the knots locations. If a single value is provided, it places the knots by means of 
  the corresponding quantile.}
  \item{order}{Spline order (cubic as default).}
  
  \item{y}{Vector of length n containing the values of the response variable.}
  \item{X}{Matrix containing the spline basis function (see details).}
  
}
\details{
  The input matrix \code{X} for \code{monspline()} can be either a B-spline basis resulting from the \code{bsb()} call, 
  or any other non-monotone spline matrix (number of rows correspond to the number of observations \eqn{n}, number of columns to 
  the number of knots). The function minimizes \eqn{((y-Xb)'(y-Xb))} over \eqn{b} that satisfy \eqn{AXb \geq 0}. \eqn{A}
  is a bidiagonal matrix of dimension \eqn{(n-1) \times n} 1's and -1's in the diagonals created internally.  
  In example section we show how to transform a B-spline into a monotone B-spline.
}

\value{
  \code{bsb} returns a list containing: 
  
  \item{X}{The B-spline basis matrix.}
  \item{knots}{Knot locations.}

  \code{monspline} returns a list containing:
  
  \item{yhat}{Fitted spline values for input basis.}
  \item{RSS}{Residual SS from the regression fit.}
}


\author{Patrick Mair, Jan de Leeuw}
\seealso{\code{\link{splinefun}}}
\examples{
set.seed(123)
n <- 20
x <- sort(runif(n, -1, 1))     ## predictor
y <- x^3 + rnorm(n, 0, 0.1)    ## response

knots <- 4                     ## knots at quartiles
bsplBasis <- bsb(x, knots)     ## B-spline basis

fit <- monspline(y, bsplBasis$X) ## monotone B-spline
fit$RSS

plot(x, y, pch = 19, main = "Monotone B-Spline")
lines(x, fit$yhat, lwd = 2, col = "coral4")
abline(v = bsplBasis$knots, col = "lightgray", lty = 2)

}
\keyword{models}
