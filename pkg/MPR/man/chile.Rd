\name{chile}
\alias{chile}
\docType{data}
\title{
Chile dataset.
}
\description{
This dataset is a modified version of the dataset used in Wright and London (2009), originally taken from pepperjoe.com. The chile length is categorized from longest to shortest. 
}
\usage{
data("chile")
}

\format{
 A data frame with 85 chiles and the following 3 variables.
  \describe{
    \item{\code{name}}{Chile name.}
    \item{\code{length}}{Chile length: ordinal (1 ... longest, 10 ... shortest)}
    \item{\code{heat}}{Chile heat scale (see details)}
}
}

\details{
Heat scale according to pepperjoe.com: 1-2 ... for sissys; 3-4 ... sort of hot; 5-6 ... fairly hot; 7-8 ... real hot; 9.5-9 ... torrid; 9.5-10 ... nuclear.
}

\source{
Wright, D. B., & London, K. (2009). Modern Regression Techniques Using R. Sage.
}


\examples{
data(chile)
str(chile)
}
\keyword{datasets}
