\name{smacof-package}
\alias{smacof-package}
\alias{smacof}
\docType{package}
\title{smacof for Multidimensional Scaling.
}
\description{This package provides various versions of multidimensional scaling (MDS) based on stress minization by means of majorization (smacof): Simple smacof on symmetric dissimilarity matrices, smacof for rectangular matrices, three-way smacof as well as extensions in terms of quadratic surfaces and simple constraints on the configuration.
}
\details{
\tabular{ll}{
Package: \tab smacof\cr
Type: \tab Package\cr
Version: \tab 0.0-1\cr
Date: \tab 2008-01-17\cr
License: \tab GPL\cr
}
The function for basic smacof on symmetric dissimilarity matrices is \code{smacofSym}.

}
\author{
Jan de Leeuw, Patrick Mair
Who wrote it

Maintainer: Jan de Leeuw <deleeuw@stat.ucla.edu>

}
\references{
de Leeuw, J. \& Mair, P. (2008). Multidimensional scaling using majorization: The R package smacof.
}

\keyword{ package }
\seealso{\code{\link{smacofSym}}

}
\examples{
data(ekman)
res <- smacofSym(ekman)

}
