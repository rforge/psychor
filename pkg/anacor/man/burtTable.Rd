\name{burtTable}
\alias{burtTable}

\title{Creates Burt Matrix}
\description{Utility function to produce a Burt matrix out of a data-frame or a profile frequency matrix
}
\usage{
burtTable(data, pf = FALSE)
}

\arguments{
  \item{data}{Data frame to be converted.}
  \item{pf}{\code{TRUE} for profile frequency matrix.}
}

\seealso{\code{\link{expandFrame}}, code{\link{mkIndiList}}}
\examples{

# artificial data 
x1 <- sample(1:6, 20, replace = TRUE) 
x2 <- sample(1:3, 20, replace = TRUE)
data <- data.frame(x1,x2)
burtTable(data)

}

\keyword{models}

