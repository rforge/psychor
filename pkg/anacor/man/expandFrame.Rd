\name{expandFrame}
\alias{expandFrame}

\title{Converts Table to Indicator Matrix}
\description{This utility function expands a matrix or data frame to an indicator supermatrix and
optionally converts this to a data frame again. By default \code{NA} becomes zero and constant rows and columns are eliminated.
}

\usage{
expandFrame(tab, clean = TRUE, zero = TRUE, returnFrame = TRUE)
}


\arguments{
  \item{tab}{Data frame of dimension n times m with frequencies. Missings are coded as \code{NA}.}
  \item{clean}{If \code{TRUE}, rows and colums with 0 margins in data frame are deleted.}
  \item{zero}{If \code{TRUE}, \code{NA}'s are replaced by 0.}
  \item{returnFrame}{If \code{TRUE}, a data frame is returned; if \code{FALSE} a matrix.}
}

\seealso{\code{\link{burtTable}}, \code{\link{mkIndiList}}}
\examples{

## Indicator matrix for tocher data
data(tocher) 
expandFrame(tocher)

}

\keyword{models}

