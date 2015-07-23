\name{plot.anacor}
\alias{plot.anacor}
\alias{plot3dstatic}
\alias{plot3dstatic.anacor}
\alias{plot3d.anacor}

\title{2-D and 3D plots for anacor solution}
\description{These functions produce various 2- and 3-dimensional plots for objects of class \code{"anacor"} 
}
\usage{
\method{plot}{anacor}(x, plot.type, plot.dim = c(1,2), legpos = "top", arrows = FALSE, conf = 0.95, 
wlines = 0, xlab, ylab, main, type, xlim, ylim, cex.axis2, ...)

\method{plot3d}{anacor}(x, plot.type, plot.dim = c(1,2,3), col.r = "RED", col.c = "BLUE", 
arrows = TRUE, xlab, ylab, zlab, main, ...)

\method{plot3dstatic}{anacor}(x, plot.type, plot.dim = c(1,2,3), col.r = "RED", col.c = "BLUE", 
arrows = TRUE, main, xlab, ylab, zlab, xlim, ylim, zlim, ...)

}

\arguments{
  \item{x}{Object of class \code{"anacor"}.}
  \item{plot.type}{Type of plot to be produced (details see below): 2-D and 3-D for \code{"jointplot"}, \code{"rowplot"}, and \code{"colplot"}; 2-D for \code{"regplot"}, \code{"graphplot"}, \code{"benzplot"}, \code{"transplot"}, and \code{"orddiag"}.}
  \item{plot.dim}{Vector of length 2 with Dimensions to be plotted. For \code{"regplot"} a single value should be provided, for  \code{"transplot"} more than two dimensions are allowed, and for \code{"benzplot"} this argument is ignored.}
  \item{legpos}{Position of the legend (for \code{"transplot"} only)}
  \item{conf}{Ellipsoid confidence level for \code{"jointplot"}, \code{"rowplot"}, and \code{"colplot"}, assuming that the ellipse where computed in \code{anacor()}. If \code{NULL}, no ellipsoids are drawn.}
  \item{arrows}{If \code{TRUE} the points in \code{"jointplot"}, \code{"rowplot"}, and \code{"colplot"} are connected with the origin.}
  \item{wlines}{For \code{"graphplot"} only: If 0, all lines are of the same thickness. For values > 0 line thickness indicates the strength of the pull}
  \item{col.r}{Color for row score plots}
  \item{col.c}{Color for column score plots}
  \item{xlab}{Label x-axis.}
  \item{ylab}{Label y-axis.}
  \item{zlab}{Label z-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{zlim}{Scale z-axis.}
  \item{main}{Plot title.}
  \item{type}{Whether points, lines or both should be plotted; for \code{"regplot"} and \code{"transplot"} only.}
  \item{cex.axis2}{For \code{"regplot"} only. The magnification to be used for the category labels in the scaled solution relative to the current setting of cex.}
  \item{\dots}{Additional graphical parameters.}
}
\details{
The following plot types are provided: \code{"jointplot"} plots row and column scores into the same device, \code{"rowplot"} and \code{"colplot"} plot the row scores and column scores, respectively, in separate devices. For these types of plots 3-dimensional versions are provided. The graph plot is an unlabeled version of the joint plot where the points are connected by lines. Options are provided (\code{wlines}) to steer the line thickness indicating the connection strength.

The regression plot (\code{"regplot"}) provides two plots. 
First, the unscaled solution is plotted. A frequency grid for the row 
categories (x-axis) and column categories (y-axis) is produced. 
The regression line is based on the category weighted means of the relative frequencies:  
the blue line on the column-wise means on the x-axis and the column category on the y-axis, 
the red line is based on the row categories on the x-axis and the row-wise means on the y-axis. 
In a second device the scaled solution is plotted. The frequency grid is determined by 
the row scores (x-axis) and the column scores(y-axis). 
Now, instead of the row/column categories, the column scores (black line y-axis) 
and the row scores (red line x-axis) are used.

The transformation plot (\code{"transplot"}) plots the row/column categories against the row/column scores. The Benzecri plot (\code{"benzplot"}) plots the observed distances against the fitted distances. It is assumed that the CA result is Benzecri scaled. The ordination diagram (\code{"orddiag"}) for CCA produces a joint plot and includes the column and row covariates based on intraset correlations. 
}

\references{De Leeuw, J. and Mair, P. (2009). Simple and Canonical Correspondence Analysis Using the R Package anacor. Journal of Statistical Software, 31(5), 1-18. \url{http://www.jstatsoft.org/v31/i05/}}
\author{Jan de Leeuw, Patrick Mair}

\seealso{\code{\link{anacor}}}

\examples{

##Regression plots using Glass data
data(glass)
res <- anacor(glass)
plot(res, plot.type = "regplot", xlab = "fathers occupation", ylab = "sons occupation")


## Benzecri Plots for bitterling data
data(bitterling)
res1 <- anacor(bitterling, ndim = 2, scaling = c("Benzecri", "Benzecri"))
res2 <- anacor(bitterling, ndim = 5, scaling = c("Benzecri", "Benzecri"))
res2
plot(res1, plot.type = "benzplot", main = "Benzecri Distances (2D)")
plot(res2, plot.type = "benzplot", main = "Benzecri Distances (5D)")

## Column score plot,transformation plot, and ordination diagram for canonical CA
data(maxwell)
res <- anacor(maxwell$table, row.covariates = maxwell$row.covariates, 
scaling = c("Goodman", "Goodman"))
res
plot(res, plot.type = "colplot", xlim = c(-1.5,1), conf = NULL)
plot(res, plot.type = "transplot", legpos = "topright")
plot(res, plot.type = "orddiag")
}

\keyword{hplot}

