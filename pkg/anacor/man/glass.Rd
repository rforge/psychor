\name{glass}
\alias{glass}
\docType{data}
\title{Glass data}
\description{Table with occupational status of fathers versus occupational status of their sons for a sample of 3497 British families.
}
\usage{data(glass)}
\format{
  Rows represent occupation of fathers, columns occupation of sons.
  \describe{
    \item{\code{PROF}}{professional and high administrative}
    \item{\code{EXEC}}{managerial and executive}
    \item{\code{HSUP}}{higher supervisory}
    \item{\code{LSUP}}{lower supervisory}
    \item{\code{SKIL}}{skilled manual and routine nonmanual}
    \item{\code{SEMI}}{semi-skilled manual}
    \item{\code{UNSK}}{unskilled manual}
  }
}

\references{
Glass, D.V. (1954). Social Mobility in Britain. Glencoe: Free Press.
}
\examples{
data(glass)
## maybe str(glass) ; plot(glass) ...
}
\keyword{datasets}
