\name{Gifi}
\alias{Gifi}
\alias{print.gifi}
\alias{summary.gifi}
\alias{scores}
\alias{scores.gifi}
\title{Gifi Methods for Categorical (nonlinear) Multivariate Analysis}
\description{
  This is the core function of the package, used by most of the wrapper functions. 
  Through this function the user can specify special models not covered by the wrapper functions. 
}
\usage{
Gifi(data, ndim = 2, rank = ndim, level = "nominal", sets = 0, active = TRUE, 
       eps = 1e-06, itermax = 1000, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data in data frame.}
  \item{ndim}{Number of dimensions to be extracted.}
  \item{rank}{Which quantification ranks. Default is number of dimensions \code{ndim}}
  \item{level}{Which quantification levels. Possible values are \code{"nominal"},  
  \code{"ordinal"}, \code{"numerical"}, and \code{"polynomial"} which can be defined as single 
  character (if all variable are of the same level) or as vector which length corresponds to the
  number of variables.} 
  \item{sets}{List of vectors of indices.}
  \item{active}{Which variables are active for computation (\code{TRUE} means all)}
  \item{eps}{Iteration precision eigenvalues.}
  \item{itermax}{Maximum number of iterations.}
  \item{verbose}{0 for no iteration output, 1 prints outer iterations, 2 prints outer and inner iterations, 3 prints outer and inner and innermost iterations.}
}
\details{
  The general idea of Gifi models is to make a joint plot in p-space of all
  objects (or individuals) and the categories of all variables. 
  
  Rank constraints require the category quantifications of the categories of
  a variable to lie in a subspace of p-space. Requiring rank equal to one for all
  variables reduces homogeneity analysis to principal component analysis
  (with optimal scaling of the variables). 
  
  Sets of variables are incorporated by using additivity restrictions on the
  category quantifications (i.e. we code the variables within a set interactively,
  but then use quantifications based on main effects only). 
  
  By combining the various types of restrictions we obtain far-reaching
  generalizations of principal component analysis, canonical analysis,
  discriminant analysis and regression analysis.
}

\value{
  Returns an object of class \code{"gifi"}.
  
  \item{datname}{Name of the data file.}
  \item{catscores}{List of category scores for each variable.}
  \item{scoremat}{Array containing the reproduced data matrix based on category scores.}
  \item{objscores}{Matrix with object scores on each dimension (scaled as z-scores).}
  \item{cat.centroids}{List of category centroids.}
  \item{ind.mat}{Dummy coded indicator matrix.}
  \item{loadings}{List of loadings of the variables on each dimension.}
  \item{low.rank}{List of lower rank quantifications.}
  \item{discrim}{Matrix with discrimination measures for each variable on each dimension.}
  \item{niter}{Number of iterations.}
  \item{eigenvalues}{Final eigenvalues.}
  \item{loss}{Value of the loss function.}
  \item{rank.vec}{Vector with ranks for each variable.}
  \item{active}{Vector with active/inactive variables.}
}
\references{
De Leeuw, J., & Mair, P. (2009). Gifi methods for optimal scaling in R: The package homals. 
Journal of Statistical Software, 31(4), 1-20, \url{http://www.jstatsoft.org/v31/i04/}. 

Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley.
}
\author{Patrick Mair, Jan de Leeuw}
\seealso{\code{\link{princals}}}
\examples{
    
  ##simple Gifi 2D-homals solution (School not active)
  res <- Gifi(galo, active = c(rep(TRUE, 4), FALSE))
}
\keyword{models}
