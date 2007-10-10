\name{plot.homals}
\alias{plot.homals}
\alias{plot3d.homals}
\alias{scatterplot3d.homals}
\title{Plot methods for homals}
\description{These methods provide various 2D and 3D plots for objects of class \code{"homals"}.
}
\usage{
plot.homals(x, plot.dim = c(1, 2), plot.type, var.subset, main, type, xlab, ylab, xlim, ylim, 
leg.pos = "topright", ...)

plot3d.homals(x, plot.dim = c(1, 2, 3), plot.type, var.subset, type, xlab, ylab, 
zlab, col, main, sphere = TRUE, bgpng = "particle.png", ax.grid = TRUE, ...)

scatterplot3d.homals(x, plot.dim = c(1, 2, 3), plot.type, var.subset, main, type, 
xlab, ylab, zlab, xlim, ylim, zlim, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"homals"}}
  \item{plot.dim}{Vector with dimensions to be plotted against.}
  \item{plot.type}{String indicating which type of plot to be produced: \code{"catplot"} 
  (also for scatterplot3d), \code{"graphplot"}, \code{"hullplot"}, \code{"labplot"} (also for scatterplot3d), 
  \code{"lossplot"}, \code{"objplot"} (also for scatterplot3d), \code{"prjplot"}, \code{"spanplot"}, \code{"starplot"}, \code{"trfplot"}, \code{"vecplot"}, 
  \code{"vorplot"}.}
  \item{var.subset}{Numeric vector for subsetting variables to be plotted. If missing, all
  variables are taken into account. Ignored for non separate variable plots.} 
  \item{main}{Plot title.}
  \item{type}{Type of points/lines to be plotted.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{zlab}{Label of z-axis.}
  \item{col}{Plotting color.}
  \item{xlim}{Limits for x-axis.}
  \item{ylim}{Limits for y-axis.}
  \item{zlim}{Limits for z-axis.}
  \item{leg.pos}{Position of the legend (for \code{"trfplot"} and \code{"spanplot"} only) 
  as provided in \code{\link{legend}}.}
  \item{sphere}{Whether 3D-sphere should be drawn.}
  \item{bgpng}{Background image in /library/rgl/textures/. If \code{NULL}, blank background.}
  \item{ax.grid}{If \code{TRUE} grid is drawn.}
  \item{\dots}{Further arguments passed to \code{\link{plot}}}
}
\details{The \code{rgl} method \code{plot3d} produces a object plot only, whereas 
\code{scatterplot3d} allows for object plots, label plots, and category quantification plots.
 
- Object plot (\code{plot.type = "objplot"}): Plots the scores of the objects (rows in data set) on two or
three dimensions.

- Category plot (\code{plot.type = "catplot"}): Plots the rank-restricted category quantifications for each
variable separately. Three-dimensional plot is available.

- Voronoi plot (\code{plot.type = "vorplot"}): Produces a category plot with Voronoi regions.

- Joint plot (\code{plot.type = "jointplot"}): Plots the object scores and the category quantifications in 1 device. 
Three-dimensional version is available.

- Graph plot (\code{plot.type = "graphplot"}): Joint plot with connections between scores/quantifications. 

- Hull plot (\code{plot.type = "hullplot"}): For each single variable the object scores are mapped onto two
dimensions and the convex hull for each response category is drawn.

- Label plot (\code{plot.type = "labplot"}): Similar to object plot, the object scores are plotted but for
each variable separately with the corresponding category labels. A three-dimensional
version is provided.

- Span plot (\code{plot.type = "spanplot"}): As label plot, it maps the object scores for each variable and
span plot connects them by the shortest path within each response category.

- Star plot (\code{plot.type = "starplot"}): Again, the object scores are mapped. In addition these points
are connected with the category centroid. Three-dimensional version is provided.

- Loss plot (\code{plot.type = "lossplot"}): Plots the rank-restricted category quantifications against the
unrestricted for each variable separately.

- Projection plot (\code{plot.type = "prjplot"}): For variables of rank 1 the category scores (two-dimensional)
are projected onto a straight line determined by the rank restricted category quantifications.

- Vector plot (\code{plot.type = "vecplot"}): For variable of rank 1 the object scores (two-dimensional) are
projected onto a straight line determined by the rank restricted category quantifications.

- Transformation plot (\code{plot.type = "trfplot"}): Plots the original (categorical) scale against the
transformed (metric) scale on each dimension over the categories of each variable separately.
}

\seealso{\code{\link{homals}}}
\examples{

##Graphplot for Hartigan solution: 3 dimensions extracted, dimension 1 plotted
##against dimension 3.
data(hartigan)
res <- homals(hartigan, ndim = 3)
plot(res, plot.dim = c(1,3), plot.type = "graphplot")

##3D star plot for indentation and bottom
plot3d(res, plot.type = "starplot", var.subset = c(2,4))

##Static 3D joint plot for thread and bottom
scatterplot3d(res, plot.type = "jointplot", var.subset = c(1,4))
}

\keyword{ hplot }
