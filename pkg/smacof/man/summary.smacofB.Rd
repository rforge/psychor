\name{summary.smacofB}
\alias{summary.smacofB}
\alias{print.smacofB}
\alias{summary.smacofR}
\alias{print.smacofR}
\alias{summary.smacofID}
\alias{print.smacofID}


\title{S3 methods for smacof}
\description{S3 methods for objects of class \code{smacofB}, \code{smacofR}, and \code{smacofID}.
}
\usage{
\method{summary}{smacofB}(object, ...)
\method{print}{smacofB}(x, ...)
\method{summary}{smacofR}(object, ...)
\method{print}{smacofR}(x, ...)
\method{summary}{smacofID}(object, ...)
\method{print}{smacofID}(x, ...)
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

