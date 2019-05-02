\name{vmu}
\alias{vmu}
\alias{print.vmu}
\alias{plot.vmu}

\title{Vector Model of Unfolding}
\description{
Computes the metric vector model of unfolding (VMU) on rectangular input data (preferences, ratings) with the individuals (rows) represented as vectors in the biplot. 
}
\usage{
vmu(delta, ndim = 2, center = TRUE, scale = FALSE)

\method{plot}{vmu}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Data frame or matrix of preferences, ratings, dissimilarities}
  \item{ndim}{Number of dimensions}
  \item{center}{If \code{TRUE} input data are centered row-wise.}
  \item{scale}{If \code{TRUE} input data are scaled row-wise.}
  \item{x}{Object of class \code{"vmu"}.}
  \item{...}{Additional arguments passed to \code{biplot} in \code{stats}.}
}

\value{
  \item{conf.row}{Row coordinates}
  \item{conf.col}{Column coordinates}
  \item{VAF}{variance accounted for}
}

\references{
Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.

Borg, I., Groenen, P. J. F., & Mair, P. (2018). Applied Multidimensional Scaling and Unfolding (2nd ed.). Springer.

Tucker, L. R. (1960). Intra-individual and inter-individual multidimensionality. In H. Gulliksen & S. Messick (Eds.), Psychological scaling: Theory and applications (pp. 155-167). Wiley.
}

\author{Ingwer Borg and Patrick Mair}

\seealso{\code{\link{biplot}}, \code{\link{unfolding}}}

\examples{
## VMU on portrait value questionnaire ratings
fit_vmu <- vmu(PVQ40agg)         ## fit 2D VMU
fit_vmu
plot(fit_vmu, cex = c(1, 0.7))   ## call biplot from stats
}

\keyword{multivariate}
