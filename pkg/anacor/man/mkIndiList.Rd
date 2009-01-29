\name{mkIndiList}
\alias{mkIndiList}

\title{Converts Table to Indicator Matrix}
\description{This function takes a data frame, a vector of types, a list of knot vectors, and a vector
of orders. It returns a list of codings for the variables, i.e., crisp indicator, numerical version,
or fuzzy indicator. 

}
\usage{
mkIndiList(data, type = rep("C",dim(data)[2]), knots = repList(0,length(type)),
ord = rep(1,length(type)))
}

\arguments{
  \item{data}{Data frame to be converted.}
  \item{type}{If \code{"C"}, a crisp indicator is returned; if \code{"A"}, a numerical version; if \code{"F"}, the b-spline basis as a fuzzy indicator.} 
  \item{knots}{Knot sequence for type-C coding.}
  \item{ord}{Variable order for type-C coding.}
}


\seealso{\code{\link{expandFrame}}, \code{\link{burtTable}}}
\examples{

# artificial data 
x1 <- sample(1:6, 20, replace = TRUE) 
x2 <- sample(1:3, 20, replace = TRUE)
data <- data.frame(x1,x2)
mkIndiList(data)

}

\keyword{models}

