\name{vss}
\alias{vss}
\alias{vss.princals}
\alias{print.princalsvss}
\alias{summary.princalsvss}
\alias{print.summary.princalsvss}
\alias{plot.princalsvss}
\title{Very Simple Structure Criterion}
\description{
  Applies the Very Simple Structure (VSS) criterion and the Velicer MAP to support the user in finding a the number of components 
  needed in a PRINCALS solution.
}
\usage{
vss(x, ndim = NULL)

\method{plot}{princalsvss}(x, main = "Very Simple Structure", 
                           xlab = "Number of Components", 
                           ylab = "VSS", ylim = c(0,1), leg.pos = "bottomleft", 
                           leg.lab = c("VSS1", "VSS2"), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"princals"}.}
  \item{ndim}{Number of dimensions. If \code{NULL}, the dimensionality of the PRINCALS fit is used. If an integer values, PRINCALS is 
  re-fitted with the corresponding number of dimensions.}
  
  \item{main}{Plot title.}
  \item{xlab}{Label x-axis.}
  \item{ylab}{Label y-axis.}
  \item{ylim}{Limits y-axis.}
  \item{leg.pos}{Legend position.}
  \item{leg.lab}{Legend labels.}
  \item{...}{Further plot arguments.}
}
\details{
  Internally it uses the \code{vss()} function from the \code{psych} package with the proper PCA settings. 
  If a rotation was used in the PRINCALS fit, the same rotation is used for VSS. 
  
}

\value{
  Returns an object of class \code{"princalsvss"}. 
  
  \item{VSS1}{VSS complexity 1.}
  \item{VSS2}{VSS complexity 2.}
  \item{MAP}{Map criterion.}
  \item{sqresid}{Squared residual correlations.}
  \item{fit}{Component fit for the complete model.}
}
\references{
Revelle, W. and Rocklin, T. (1979). Very Simple Structure: an alternative procedure for estimating 
the optimal number of interpretable factors, Multivariate Behavioral Research, 14, 403-414. 

Velicer, W. (1976) Determining the number of components from the matrix of partial correlations. Psychometrika, 41, 321-327.
}
\author{Patrick Mair, Jan de Leeuw}
\seealso{\code{\link{princals}}, \code{\link[psych]{vss}}}
\examples{
## extract subset of ABC data
ABC6 <- ABC[, 6:11]

## ordinal PCA, 2 dimensions
fit <- princals(ABC6)

## VSS for 5 dimensional solution
fitvss <- vss(fit, ndim = 5)
fitvss
summary(fitvss)
plot(fitvss)
}
\keyword{models}
