\name{gravity}
\alias{gravity}

\title{Gravity dissimilarities}
\description{Computes the dissimilarities using an (exponential) gravity model based on co-occurences
}

\usage{
gravity(X, lambda = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A numeric matrix}
  \item{lambda}{Tuning parameter}
}

\details{The first step in this function is to compute the co-occurences. Based on the 
binarized data matrix \eqn{Y} we compute \eqn{Y'Y} which leads to the co-occurence matrix.
We then use the gravity model to compute the gravity dissimilarities. 
In order to give more (or less) structure to the MDS solution, the tuning parameter lambda can 
be increased (or decresed with lower bound 0). In addition, a weight matrix is created 
that sets cells with no co-occurences to 0, i.e. they are blanked out in a subsequent 
smacof fit.
}

\value{
  \item{gravdiss}{Gravity dissimilarities}
  \item{weightmat}{Weight matrix for subsequent smacof computation}
  \item{co.occ}{Matrix with co-occurences}
}
\references{Mair, P., Rusch, T. & Hornik, K. (2014). Gravity MDS on sparse document-term matrices: 
Republican values. Research Report Series, Department of Statistics and Mathematics, WU
Vienna, Report 127.
}
\author{Patrick Mair}

\seealso{\code{\link{smacofSym}}}

\examples{
data(GOPdtm)
gravD <- gravity(GOPdtm, lambda = 2)
res <- smacofSym(gravD$gravdiss, weightmat = gravD$weightmat)
plot(res)
}

\keyword{models}
