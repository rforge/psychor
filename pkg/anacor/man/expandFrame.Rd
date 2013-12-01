\name{expandFrame}
\alias{expandFrame}

\title{Expand Matrix}
\description{This utility function expands a matrix or data frame to an indicator supermatrix and
optionally converts this to a data frame again. By default \code{NA} becomes zero and constant rows and columns are eliminated.
}

\usage{
expandFrame(tab, clean = TRUE, zero = TRUE, returnFrame = TRUE)
}


\arguments{
  \item{tab}{Data frame (factors). Missings are coded as \code{NA}.}
  \item{clean}{If \code{TRUE}, rows and colums with 0 margins in data frame are deleted.}
  \item{zero}{If \code{TRUE}, \code{NA}'s are replaced by 0.}
  \item{returnFrame}{If \code{TRUE}, a data frame is returned; if \code{FALSE} a matrix.}
}

\seealso{\code{\link{burtTable}}, \code{\link{mkIndiList}}}
\examples{

## sleeping bags
data(sleeping)
sleeping_cat <- sleeping
temp_cat <- cut(sleeping$Temperature, c(-20, -1, 7), labels = c("warm", "cold")) 
sleeping_cat$Temperature <- temp_cat
weight_cat <- cut(sleeping$Weight, c(700, 1100, 2200), labels = c("light", "heavy")) 
sleeping_cat$Weight <- weight_cat
price_cat <- cut(sleeping$Price, c(100, 250, 350, 700), 
labels = c("cheap", "medium", "expensive"))  
sleeping_cat$Price <- price_cat
sleeping_cat
expandFrame(sleeping_cat)

}

\keyword{models}

