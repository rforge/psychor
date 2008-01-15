\name{anacor-package}
\alias{anacor-package}
\docType{package}
\title{Simple and Canonical Correspondence Analysis
}
\description{
This package performs simple correspondence analysis (CA) on a two-way frequency table (with missings) by means of SVD. Different scaling methods (standard, centroid, Benzecri, Goodman) as well as various plots including confidence ellipsoids are provided. }
\details{
\tabular{ll}{
Package: \tab anacor\cr
Type: \tab Package\cr
Version: \tab 0.9-0\cr
Date: \tab 2007-12-19\cr
License: \tab GPL\cr
}
The main function of the package is \code{anacor()} which performs simple CA. The input dataset \code{tab} is a (incomplete) frequency table. 
Various scaling options in terms of standard, centroid, Benzecri, and Goodman scaling are provided. They can be specified using the \code{scaling} argument.
If covariates for row/column scores are provided, canonical CA will be computed. 
Utility functions to create various types of indicator matrices are implemented.
}
\author{Jan de Leeuw, Patrick Mair

Maintainer: Jan de Leeuw <deleeuw@stat.ucla.edu>
}
\references{
de Leeuw, J. and Mair, P. (2007). Simple and Canonical Correspondence Analysis Using the R Package anacor. 
}

\keyword{ package }
\seealso{\code{\link{anacor}}, \code{\link[<pkg>:<pkg>-homals]{<pkg>}}, \code{\link[<pkg>:<pkg>-ca]{<pkg>}}, \code{\link[<pkg>:<pkg>-ade4]{<pkg>}}
}
\examples{
data(tocher)
res <- anacor(tocher)
res
summary(res)
}
