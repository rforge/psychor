\name{plot.princals}
\alias{plot.princals}

\title{Plot Methods for Princals}
\description{Various plotting options for objects of class \code{"princals"}.
}
\usage{
\method{plot}{princals}(x, plot.type = "biplot", plot.dim = c(1, 2), var.subset = "all", 
col.scores = "gray", col.loadings = "black", col.lines = "black", cex.scores = 0.8, 
cex.loadings = 0.8, labels.scores = FALSE, labels.loadings = TRUE, asp = 1, main, xlab, 
ylab, xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"princals"}}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{plot.type}{String indicating which type of plot to be produced: \code{"biplot"}, 
  \code{"loadplot"}, \code{"screeplot"}, \code{"transplot"}. See details.} 
  \item{var.subset}{Numeric vector for subsetting variables to be plotted. If missing, all variables are taken into account. Ignored for non separate variable plots.} 
  \item{col.scores}{Color object scores.} 
  \item{col.loadings}{Color loading vectors.}
  \item{col.lines}{Line colors for transformation and scree plot.}
  \item{cex.scores}{Size object scores.}
  \item{cex.loadings}{Size loading labels.}
  \item{asp}{Aspect ratio}
  \item{labels.scores}{Whether object score labels should be plotted.} 
  \item{labels.loadings}{Whether variable names should be plotted.}
  \item{main}{Plot title.}
  \item{type}{Type of points/lines to be plotted.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Limits for x-axis.}
  \item{ylim}{Limits for y-axis.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}}
}
}
\details{Plot description:
 
- Biplot (\code{plot.type = "biplot"}): Biplot with the object scores and the loadings.

- Transformation plot (\code{plot.type = "transplot"}): Plots the original scale against the
transformed scale (first dimension) for each variable separately.

- Loadings plot (\code{plot.type = "loadplot"}): Plots the loadings of the variables and connects them with the origin. 

- Scree plot (\code{plot.type = "screeplot"}): Produces a scree plot based on the eigenvalues. 
}

\examples{
ABC6 <- ABC[,6:11]
fitord <- princals(ABC6, ndim = 3)  ## 3D princals solution

## transformation plot
plot(fitord, plot.type = "transplot", var.subset = 1:4, lwd = 2, col.lines = "coral")

## loadings plot
plot(fitord, "loadplot", main = "Loadings Plot ABC Data")

## biplot
plot(fitord, "biplot", labels.scores = TRUE, main = "Biplot ABC Data", cex.scores = 0.5)  
plot(fitord, "biplot", main = "Biplot ABC Data")
plot(fitord, "biplot", main = "Biplot ABC Data", plot.dim = c(2,3))

## scree plot
plot(fitord, "screeplot")  
}

\keyword{ hplot }