\name{morals}
\alias{morals}
\alias{print.morals}
\alias{predict.morals}
\alias{print.predict.morals}
\alias{summary.predict.morals}
\title{Multiple Montone Regression}
\description{
  Performs a monotone multiple regression. 
}
\usage{
morals(data, resp, level = "ordinal", active = TRUE, 
       eps = 1e-06, itermax = 1000, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data in data frame.}
  \item{resp}{Which column of the data should be considered as response variable}
  \item{level}{Which quantification levels. Possible values are \code{"nominal"},  
  \code{"ordinal"}, \code{"numerical"}, and \code{"polynomial"} which can be defined as single 
  character (if all variable are of the same level) or as vector which length corresponds to the
  number of variables.} 
  \item{active}{Which variables are active for computation (\code{TRUE} means all)}
  \item{eps}{Iteration precision eigenvalues.}
  \item{itermax}{Maximum number of iterations.}
  \item{verbose}{0 for no iteration output, 1 prints outer iterations, 2 prints outer and inner iterations, 3 prints outer and inner and innermost iterations.}
}
\details{
  This function extends simple monotone regression analysis to multiple predictors. 
  Internally the \code{gifi()} function is used with a corresponding specification 
  of the variable sets: the response goes into the first set, the predictors into the second set.
}

\value{
  Returns an object of class \code{"morals"}.
  
  \item{datname}{Name of the data file.}
  \item{catscores}{List of category scores for each variable.}
  \item{scoremat}{Array containing the reproduced data matrix based on category scores.}
  \item{objscores}{Matrix with object scores on each dimension.}
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
    
  ## MORALS solution
  res <- morals(neumann, resp = 3, verbose = 1)
  
}
\keyword{models}
