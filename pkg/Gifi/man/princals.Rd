\name{princals}
\alias{princals}
\alias{print.princals}
\alias{summary.princals}
\alias{print.summary.princals}
\alias{plot.princals}
\title{Categorical Principal Component Analysis}
\description{
  Performs categorical PCA (aka nonlinear PCA). Internally it uses \code{gifi()} with a rank-1 restriction. 
  It presents and scales eigenvalues, explained amount of variance, and loadings such that they are comparable 
  to standard PCA as e.g. fitted using \code{princomp()}. Compared to standard PCA, categorical PCA abandons the 
  assumption of metric input data as well as linearity among the observed variables.
}
\usage{
princals(data, ndim = 2, level = "ordinal", active = TRUE, eps = 1e-6, itermax = 1000, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input data as data frame.}
  \item{ndim}{Number of dimensions to be extracted.}
  \item{level}{Which quantification levels. Possible values are \code{"nominal"},  
  \code{"ordinal"}, \code{"numerical"}, and \code{"polynomial"} which can be defined as single 
  character (if all variable are of the same level) or as vector which length corresponds to the
  number of variables.} 
  \item{active}{Which variables are active for computation (\code{TRUE} means all)}
  \item{eps}{Iteration precision.}
  \item{itermax}{Maximum number of iterations.}
  \item{verbose}{0 for no iteration output, 1 prints outer iterations, 2 prints outer and inner iterations, 3 prints outer and inner and innermost iterations.}
}
\details{
  Details about the computation of PCA-comparable eigenvalues, explained amount of variance, and loadings can be 
  found in the package vignette. 
}

\value{
  Returns an object of class \code{"princals"}. The specific output values are below, the remaining ones are the 
  same as in \code{\link{Gifi}}.
  
  \item{eigenvalues}{Eigenvalues.}
  \item{vaf}{variance accounted for (component-wise in \%).}
  \item{loadings}{Loadings.}
  \item{pcscores}{Principal component scores.}
}
\references{
De Leeuw, J., & Mair, P. (2009). Gifi methods for optimal scaling in R: The package homals. 
Journal of Statistical Software, 31(4), 1-20, \url{http://www.jstatsoft.org/v31/i04/}. 

Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley.
}
\author{Patrick Mair, Jan de Leeuw}
\seealso{\code{\link{Gifi}}}
\examples{
## extract subset of ABC data
ABC6 <- ABC[, 6:11]

## --------- ordinal PCA, 1 dimension --------------------
fit2 <- princals(ABC6, ndim = 1)
fit2
summary(fit2)
fit2$loadings                         ## loadings

## comparison with standard PCA 
as.numeric.factor <- function(x) as.numeric(levels(x))[x]
ABC6m <- sapply(ABC6, as.numeric.factor)
fitpc <- princomp(ABC6m)
fitpc
fit2$eigenvalues/(fitpc$sdev^2)[1]       ## eigenvalue ratio PRINCALS/PCA

## ---- 2D solution 
fit3 <- princals(ABC6, ndim = 2)
fit3
summary(fit3)
plot(fit3)                               ## loadings plot
newdata <- scores(fit3)                  ## extract the rescaled data matrix (first dimension)
head(newdata)
## can act as input data for subsequent models

## --- Full solution
fit4 <- princals(ABC6, ndim = 6)
fit4
summary(fit4)
plot(fit4, plot.type = "screeplot")      ## scree plot

## --- Metric PCA 
fit5 <- princals(ABC6, ndim = 1, level = "numerical")
fit5
fit5$eigenvalues/(fitpc$sdev^2)[1]       ## eigenvalue ratio PRINCALS/PCA

## --- Nominal PCA
fit6 <- princals(ABC6, ndim = 1, level = "nominal")
fit6

}
\keyword{models}
