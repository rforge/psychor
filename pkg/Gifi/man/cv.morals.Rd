\name{cv.morals}
\alias{cv.morals}
\alias{cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-one-out cross-validation for Morals fit
}
\description{
Performs a LOOCV on Morals objects and returns the prediction RMSE.
}
\usage{
cv(object, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"morals"}}
  \item{verbose}{If \code{TRUE}, iteration is printed out}
  \item{...}{Additional arguments (ignored here)}
}
%- maybe also 'usage' for other objects documented here.


\examples{
## mimic linear regression
xknots <- knotsGifi(neumann[,1:2], "E")  
yknots <- knotsGifi(neumann[,3], "E")
fitlin <-  morals(neumann[,1:2], neumann[,3], xknots = xknots, yknots = yknots, xdegrees = 1, ydegrees = 1,
   xordinal = FALSE, yordinal = FALSE)
cv(fitlin)  ## RSME prediction

}
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")
