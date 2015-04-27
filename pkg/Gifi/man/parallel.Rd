\name{parallel}
\alias{parallel}
\alias{parallel.princals}
\alias{print.princalspa}
\alias{plot.princalspa}
\title{Parallel Analysis for PRINCALS}
\description{
  Compares the eigenvalue structure of a PRINCALS solution with either random ("parallel") matrices or 
  with permutated matrices. 
}
\usage{
parallel(x, N = 1000, method = "permutation")

\method{plot}{princalspa}(x, main = "Parallel Analysis", xlab = "Number of Components", 
                          ylab = "Eigenvalues", leg.pos = "topright", leg.lab = c("Fitted", "Simulated"), ylim, ...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"princals"}.}
  \item{N}{Number of simulated matrices.}
  \item{method}{Either \code{"permutation"} or \code{"random"} (see details).}
  
  \item{main}{Plot title.}
  \item{xlab}{Label x-axis.}
  \item{ylab}{Label y-axis.}
  \item{ylim}{Limits y-axis.}
  \item{leg.pos}{Legend position.}
  \item{leg.lab}{Legend labels.}
  \item{...}{Further plot arguments.}
}
\details{
  If \code{method = "random"}, \code{N} random matrices are generated of the same size as the data matrix. 
  The values are drawn from a standard normal distribution. If \code{method = "permutation"}, \code{N} 
  column-wise permutations of the PRINCALS score matrix are performed. In either case, the main output 
  are the mean of the eigenvalues across \code{N} simulated matrices.
  
  The criterion for evaluating the number of components is to check how many PRINCALS eigenvalues are 
  larger than the mean eigenvalues of the simulated matrices. 
}

\value{
  Returns an object of class \code{"princalspa"}. 
  
  \item{pceigen}{Eigenvalues of the fitted PRINCALS solution.}
  \item{simeigen}{Average eigenvalues of the simulated matrices.}
  \item{simsd}{Standard deviation of the simulated eigenvalues.}
  \item{simci}{95\% confidence intervals of simulated eigenvalues.}
  \item{ncomp}{Suggested number of components.}
  
}
\references{
Horn, J. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179-185.
}

\author{Patrick Mair, Jan de Leeuw}
\seealso{\code{\link{princals}}, \code{\link{vss}}}
\examples{
## extract subset of ABC data
ABC6 <- ABC[, 6:11]

## ordinal PCA, full dimensions
fit <- princals(ABC6, ndim = 6)

fitpar <- parallel(fit) 
fitpar
summary(fitpar)
plot(fitpar)
}
\keyword{models}
