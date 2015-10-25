\name{anacor}
\alias{anacor}
\alias{print.anacor}

\title{Simple and Canonical Correspondence Analysis}
\description{This function performs simple and canonical CA for possibly incomplete tables based on SVD. 
Different scaling methods for row and column scores are provided. 
}
\usage{
anacor(tab, ndim = 2, row.covariates, col.covariates, scaling = c("standard","standard"), 
ellipse = TRUE, eps = 1e-06)

\method{print}{anacor}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{Data frame of dimension n times m with frequencies. Missings are coded as \code{NA}.}
  \item{ndim}{Number of dimensions.}
  \item{row.covariates}{Matrix with n rows containing covariates for the row scores.}
  \item{col.covariates}{Matrix with m rows containing covariates for the column scores.}
  \item{scaling}{A vector with two elements. The first one corresponds to the method for row scaling, the second one for column scaling. Available scaling methods are \code{"standard"}, \code{"Benzecri"}, \code{"Goodman"}.}
  \item{ellipse}{If \code{TRUE}, confidence ellipses are computed.}
  \item{eps}{Convergence criterion for reconstitution algorithm.}
  \item{x}{Object of class \code{"anacor"} in \code{print.anacor}.}
  \item{\dots}{Additional arguments ignored.}
}
\details{
  Missing values in \code{tab} are imputed using the reconstitution algorithm. Setting \code{scaling} to \code{"standard"} leads to standard coordinates. Principal coordinates can be computed by means of Benzecri decomposition. Furthermore, scores can be scaled around their centroid. 
  Goodman scaling is based on Fisher-Maung decomposition.  
  
  For large datasets it is suggested to set \code{ellipse = FALSE}. If \code{ellipse = TRUE}, make sure
  that there are no rows and columns that have full 0 entries. 
}
\value{
  \item{row.scores}{Scaled row scores.}
  \item{col.scores}{Scaled column scores.}
  \item{ndim}{Number of dimensions extracted.}
  \item{chisq}{Total chi-square value.}
  \item{chisq.decomp}{Chi-square decomposition across dimensions with p-values.}
  \item{singular.values}{Singular values without trivial solution.}
  \item{se.singular.values}{Standard errors for the singular values.}
  \item{stestmat}{Test results for singular values (only for \code{ellipse = TRUE}).}
  \item{left.singvec}{Left singular vectors without trivial solution.}
  \item{right.singvec}{Right singular vectors without trivial solution.}
  \item{eigen.values}{Eigenvalues for the fitted dimensionality.}
  \item{eigenvall}{Full vector of eigenvalues (principal inertias).}
  \item{datname}{Name of the dataset.}
  \item{tab}{Table with imputed frequencies in case of missings.}
  \item{row.covariates}{Matrix with row covariates.}
  \item{col.covariates}{Matrix with column covariates.}
  \item{scaling}{Scaling Method.}
  \item{bdmat}{List of matrices with observed and fitted Benzecri distances for rows and columns.}
  \item{rmse}{Root mean squared error of Bezencri distances (rows and columns).}
  \item{row.acov}{Covariance matrix for row scores.}
  \item{col.acov}{Covariance matrix for column scores.}
  \item{cancoef}{List containing canonical coefficients (CCA only).}
  \item{sitescores}{List containing the site scores (CCA only).}
  \item{isetcor}{List containing the intraset correlations (CCA only).}
}
\references{De Leeuw, J. and Mair, P. (2009). Simple and Canonical Correspondence Analysis Using the R Package anacor. Journal of Statistical Software, 31(5), 1-18. \url{http://www.jstatsoft.org/v31/i05/}
}
\author{Jan de Leeuw, Patrick Mair}

\seealso{\code{\link{plot.anacor}}}
\examples{

## simple CA on Tocher data, symmetric standard coordinates
data(tocher)
res <- anacor(tocher)
res

## simple CA on Tocher data, asymmetric coordinates
res <- anacor(tocher, scaling = c("standard", "Benzecri"))
res

## 2- and 5-dimensional solutions for bitterling data, Benzecri scaling
data(bitterling)
res1 <- anacor(bitterling, ndim = 2, scaling = c("Benzecri", "Benzecri"))
res2 <- anacor(bitterling, ndim = 5, scaling = c("Benzecri", "Benzecri"))
res1
res2

## Canonical CA on Maxwell data, Goodman scaling
data(maxwell)
res <- anacor(maxwell$table, row.covariates = maxwell$row.covariates, 
scaling = c("Goodman", "Goodman"))
res

}

\keyword{models}

