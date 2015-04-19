\name{Gifi-package}
\alias{Gifi-package}
\docType{package}
\title{Gifi Methods for Categorical Multivariate Analysis}

\description{Gifi implements methods for multivariate analysis of categorical data. It replaces the homals package. 
The core function is called gifi and several wrapper functions are available for computing 
the following Gifi models: multiple correspondence analysis aka homogeneity analysis (HOMALS), categorical PCA (PRINCALS), multiple montone regression (MORALS).  
}
\details{
\tabular{ll}{
Package: \tab Gifi\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2015-04-19\cr
License: \tab GPL-2\cr
URL: \tab \url{http://r-forge.r-project.org/projects/psychor/}\cr
}
The basic function of the Gifi package is \code{gifi} which returns an object of class \code{"gifi"}. It allows for 
grouping variables into sets, constraining ranks on category quantifications and levels.
}
\author{
Patrick Mair, Jan de Leeuw 

Maintainer: Patrick Mair <mair@fas.harvard.edu>
}
\references{
De Leeuw, J., & Mair, P. (2009). Gifi methods for optimal scaling in R: The package homals. 
Journal of Statistical Software, 31(4), 1-20, \url{http://www.jstatsoft.org/v31/i04/}. 

Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley.
}

\keyword{ package }
\seealso{\code{\link{gifi}}, \code{\link{princals}}
}