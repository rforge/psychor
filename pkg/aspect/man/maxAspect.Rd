\name{maxAspect}
\alias{maxAspect}
\alias{print.aspect}
\alias{summary.aspect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scaling by Maximizing Correlational Aspects}
\description{This function performs optimal scaling by maximizing a certain aspect of the correlation matrix.
}
\usage{
maxAspect(data, aspect = aspectSum, itmax = 100, eps = 1e-06, extra = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame or matrix}
  \item{aspect}{Function on the correlation matrix (see details)}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{extra}{Extra parameter for aspect}
}
\details{We provide various pre-specified aspects: \code{"aspectAbs"} takes the sum of the absolute values of the correlations, \code{"aspectSum"} the sum of the correlations, \code{"aspectDeterminant"} the determinant of the correlation matrix, \code{"aspectEigen"} the sum of the first p eigenvalues (principal component analysis), \code{"aspectSMC"} the squared multiple correlations (multiple regression), \code{"aspectSumSMC"} uses the sum of p SMC's (path analysis). optionally, the user can specify his own aspects. 
}
\value{
  \item{loss}{Final value of the loss function}
  \item{catscores}{Resulting category scores (after optimal scaling)}
  \item{cormat}{Correlation matrix based on the scores}
  \item{eigencor}{Eigenvalues of the correlation matrix}
  \item{indmat}{Indicator matrix (dummy coded)}
  \item{scoremat}{Transformed data matrix (i.e with category scores resulting from optimal scaling)}
  \item{burtmat}{Burt matrix}
  \item{niter}{Number of iterations}
}
\references{
Mair, P., & de Leeuw, J. (2008). Scaling variables by optimizing correlational and non-correlational aspects in R. Journal of Statistical Software, forthcoming.

de Leeuw, J. (1988). Multivariate analysis with optimal scaling. In S. Das Gupta and J.K. Ghosh, Proceedings of the International Conference on Advances in Multivariate Statistical Analysis, pp. 127-160. Calcutta: Indian Statistical Institute.
}
\author{Jan de Leeuw, Patrick Mair}

\seealso{\code{\link{lineals}}}
\examples{
data(galo)
res.eig <- maxAspect(galo[,1:4], aspect = "aspectEigen")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

