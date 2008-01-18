\name{summary.smacofB}
\alias{summary.smacofB}
\alias{print.smacofB}

\title{S3 methods for smacof}
\description{S3 methods for objects of class \code{smacofB}
}
\usage{
\method{summary}{smacofB}(object, ...)
\method{print}{smacofB}(x, ...)
}
\arguments{
  \item{object}{Object of class \code{smacofB}}
  \item{x}{Object of class \code{smacofB}}
  \item{\dots}{Ignored}
}
\details{
to be added  
}

\examples{
data(ekman)
res <- smacofSym(ekman)
res
summary(res)

}
\keyword{models}

