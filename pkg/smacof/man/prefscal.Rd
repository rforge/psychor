\name{prefscal}
\alias{prefscal}

\title{Nonmetric unfolding}
\description{Variant of smacof for rectangular matrices (typically ratings, preferences) that allows for nonmetric transformations. Also known as nonmetric unfolding.
}
\usage{
prefscal(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
         conditionality = c("matrix", "row"), lambda = 0.5, omega = 1, 
         circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
         ties = c("primary","secondary", "tertiary"), verbose = FALSE, itmax = 1000,  
         reg = 1e-6, eps = 1e-6, spline.degree = 2, spline.intKnots = 2)
}

\arguments{
  \item{delta}{Data frame or matrix of preferences, ratings, dissimilarities.}
  \item{ndim}{Number of dimensions}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{conditionality}{A single transformations are applied for the entire matrix \code{"matrix"} or for each row separately \code{"row"}}
  \item{lambda}{Penalty strength balancing the loss contribution of stress and the penalty}
  \item{omega}{Penalty width determines for what values of the variation coefficient the penalty should become active.}
  \item{circle}{If \code{"column"}, the column configurations are restricted to be on a circle, if \code{"row"}, row configurations are on a circle, if \code{"none"}, there are no restrictions on row and column configurations}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{ties}{Tie specification for non-metric MDS}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{itmax}{Maximum number of iterations}
  \item{reg}{Regularization factor, prevents distances from being 0}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}
\details{Creates an object of class \code{smacofR}. 
}

\value{
  \item{obsdiss}{Observed dissimilarities, corresponds to delta}
  \item{confdist}{Configuration dissimilarities}
  \item{conf.row}{Matrix of final row configurations}
  \item{conf.col}{Matrix of final column configurations}
  \item{stress}{Final, normalized stress value}
  \item{spp.row}{Stress per point, rows}
  \item{spp.col}{Stress per point, columns}
  \item{congvec}{Vector of congruency coefficients}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nind}{Number of individuals (rows)}
  \item{nobj}{Number of objects (columns)}
}
\references{de Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 

Busing, F. M. T. A., Groenen, P. J. F., & Heiser, W. J. (2005). Avoiding degeneracy in multidimensional unfolding by penalizing on the coefficient of variation. Psychometrika, 70, 71-98.
}

\author{Patrick Groenen, Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{plot.smacof}}, \code{\link{smacofConstraint}}, \code{\link{smacofSym}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}}

\examples{
res <- prefscal(breakfast)
res

## various configuration plots
plot(res)
plot(res, type = "p", pch = 25)
plot(res, type = "p", pch = 25, col.columns = 3, 
     label.conf.columns = list(label = TRUE, pos = 3, col = 3), 
     col.rows = 8, label.conf.rows = list(label = TRUE, pos = 3, col = 8))

## Shepard plot
plot(res, "Shepard")

## Stress decomposition chart
plot(res, "stressplot")

}

\keyword{models}
