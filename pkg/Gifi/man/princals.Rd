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
princals(data, ndim = 2, ordinal = TRUE, ties = "s", knots = knotsD(data), 
degrees = 2, copies = 1, missing = "m", normobj.z = TRUE, active = TRUE,
itmax = 1000, eps = 1e-06, seed = 123, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input data frame: n observations, m variables}
  \item{ndim}{Number of dimensions to be computed}
  \item{ordinal}{Whether variables should be considered as ordinal or not. Alternatively, one can specify a boolean vector of length m denoting which variables should be ordinally restricted or not.}
  \item{ties}{How ties should be handled: primary (\code{"p"}), secondary (\code{"s"}), or tertiary (\code{"t"})}
  \item{knots}{Knots specification for spline transformation (see details)}
  \item{degrees}{Spline degrees. If different degrees should be used across variables, a vector of length m can be specified.}
  \item{copies}{Number of copies for each variables (also as vector of length m)}
  \item{missing}{How missing values should be handled: multiple (\code{"m"}), single (\code{"s"}), or average (\code{"a"})}
  \item{active}{Which variables should be active or inactive (also as vector of length m)}
  \item{normobj.z}{If \code{TRUE}, object scores are z-scores, if \code{FALSE}, they are restriction to SS of 1.}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{seed}{Random number seed (for starting solution of the object scores)}
  \item{verbose}{Iteration printout}
}
\details{
Details will follow.
}
\value{
\item{transform}{Optimally transformed scores}
\item{rhat}{Induced correlation matrix}
\item{evals}{Eigenvalues of induced correlation matrix}
\item{objectscores}{Object scores (rows)}
\item{scoremat}{Optimally scaled data matrix (first dimension)}
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
ABC6 <- ABC[,6:11]

## ordinal PCA
fitord <- princals(ABC6)  ## ordinal PCA
fitord
summary(fitord)

## linear restrictions (mimics standard PCA)
abc_knots <- knotsE(ABC6)     ## 0 interior knots
fitlin <- princals(ABC6, knots = abc_knots, degrees = 1)  
fitlin 
fitlin$evals

## compare with standard PCA
ABCnum <- makeNumeric(ABC6)
fitpca <- prcomp(ABCnum, scale = TRUE)
fitpca$sdev^2   
}
\keyword{multivariate}% use one of  RShowDoc("KEYWORDS")
