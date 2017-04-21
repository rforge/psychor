\name{WilPat}
\alias{WilPat}
\docType{data}
\title{
Wilson-Patterson Conservatism Scale
}
\description{
This dataset contains a modified version of the classical Wilson-Patterson conservatism scale. Each item has the following response categories: 0 ... disapprove, 1 ... approve, 2 ... don't know. 
}

\usage{
data("WilPat")
}


\format{
The first 15 items are conservative items, the remaining ones are liberal. There are 804 persons in the sample. In addition there are the following covariates:
\describe{
	\item{\code{Country}}{Participant's country.}
  \item{\code{LibCons}}{Self-reported liberalism/conservatism.}
  \item{\code{LeftRight}}{Self-reported left/right identification.}
  \item{\code{Gender}}{Gender.}
  \item{\code{Age}}{Age.}
}
}

\note{
Thanks to Benedek Kurdi and Levente Littvay for providing this dataset. 
}

\examples{
data("WilPat")
str(WilPat)
}
\keyword{datasets}
