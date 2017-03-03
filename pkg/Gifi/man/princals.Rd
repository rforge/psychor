\name{princals}
\alias{princals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Categorical principal component analysis (PRINCALS). 
}
\description{
Fits a categorical PCA. 
}
\usage{
princals(data, ndim = 2, ordinal = TRUE, ties = "s", knots = knotsQ(data), degrees = 2, copies = 1, missing = "m", active = TRUE, itmax = 1000, eps = 1e-06, seed = 123, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input data frame: n observations, m variables}
  \item{ndim}{Number of dimensions to be computed}
  \item{ordinal}{Whether variables should be considered as ordinal or not (can also be a boolean vector of length m}
  \item{ties}{How ties should be handled: primary, secondary, or tertiary}
  \item{knots}{Knots specification for spline transformation (see details)}
  \item{degrees}{Spline degrees}
  \item{copies}{Number of copies for each variables (can also be a vector of length m)}
  \item{missing}{How missing values should be handled}
  \item{active}{Which variables should be active or inactive (vector of length m)}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{seed}{Random number seed}
  \item{verbose}{Iteration prinout (boolean)}
}
\details{
Details will follow.
}
\value{
\item{transform}{Optimally transformed scores}
\item{rhat}{Induced correlation matrix}
\item{evals}{Eigenvalues of induced correlation matrix}
\item{objectscores}{Object scores (rows)}
\item{scoremat}{Optimally scaled data matrix}
\item{quantifications}{Category quantifications}
\item{dmeasures}{Discimination matrices}
\item{lambda}{Average discrimination matrix}
\item{weights}{Component weights}
\item{loadings}{Component loadings}
\item{ntel}{Number of iterations}
\item{f}{Loss function value}
\item{data}{Original data frame}
\item{datanum}{Numerical data frame}
\item{ndim}{Number of extracted dimensions}
\item{call}{Function call}
}
\references{
Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley.
}

\examples{
fit <- princals(ABC[,6:11])
fit
}
\keyword{models}% use one of  RShowDoc("KEYWORDS")
