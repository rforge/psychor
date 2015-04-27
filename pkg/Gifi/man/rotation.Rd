\name{rotation}
\alias{rotation}
\alias{rotation.princals}
\alias{print.princalsrot}
\title{PRINCALS Rotation}
\description{
  Performs a rotation (VARIMAX or PROMAX) after a PRINCALS fit.
}
\usage{
\method{rotation}{princals}(x, method = "varimax", normalize = TRUE, eps = 1e-5, m = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"princals"}.}
  \item{method}{Rotation method: \code{"varimax"} for VARIMAX rotation, \code{"promax"} for PROMAX rotation.}
  \item{normalize}{Logical. Should Kaiser normalization be performed?}
  \item{eps}{The tolerance for stopping: the relative change in the sum of singular values.}
  \item{m}{The power used the target for promax. Values of 2 to 4 are recommended.}
}
\details{
  Note that VARIMAX and PROXMAX rotations have been proposed within a factor-analytic context. They use factor loadings
  as input. PRINCALS loadings (or PCA loadings) in general are not factor loadings, they need to be re-scaled by their 
  respective eigenvalues before performing the rotation. The \code{rotation()} function does that internally and then uses 
  the \code{varimax()} and \code{promax()} functions from the \code{stats} package. After rotation it scales back 
  the loadings to the PRINCALS metric. As a result, the loadings returned from \code{rotation()} are on the same scale 
  as the PRINCALS loadings; the corresponding principal scores as well. This makes the rotated solution comparable 
  to the PRINCALS solution.
}

\value{
  Returns an object of class \code{"princalsrot"}. 
  
  \item{loadings}{Loadings after rotation.}
  \item{pcscores}{Principal component scores after rotation.}
  \item{rotmat}{Rotation matrix.}
}
\references{
Hendrickson, A. E., & White, P. O. (1964) Promax: a quick method for rotation to orthogonal oblique structure. 
British Journal of Statistical Psychology, 17, 65-70.

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. 
Psychometrika 23, 187-200.
}
\author{Patrick Mair, Jan de Leeuw}
\seealso{\code{\link{varimax}}, \code{\link{promax}}}
\examples{
## extract subset of ABC data
ABC6 <- ABC[, 6:11]

## ordinal PCA, 2 dimensions
fit <- princals(ABC6)

## varimax rotation
rot1 <- rotation(fit)
rot1

## promax rotation
rot2 <- rotation(fit, method = "promax")
rot2
}
\keyword{models}
