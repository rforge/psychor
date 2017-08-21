\name{FamilyIQ}
\alias{FamilyIQ}
\docType{data}
\title{
Family Intelligence
}
\description{
Dataset from Hox (2010) containing six intelligence measures. Children are nested within families. 
}

\usage{
data("FamilyIQ")
}

\format{
A data frame with 399 children, nested within 60 families: 
  \describe{
	 \item{\code{family}}{Family ID.}
   \item{\code{child}}{Child ID.}
	 \item{\code{wordlist}}{Word list intelligence measure.}
	 \item{\code{cards}}{Cards intelligence measure.}
	 \item{\code{matrices}}{Matrices intelligence measure.}
	 \item{\code{figures}}{Figures intelligence measure.}
	 \item{\code{animals}}{Animals intelligence measure.}
	 \item{\code{occupation}}{Occupation intelligence measure.}
}
}

\source{
Hox, J. J. (2010). Multilevel analysis: Techniques and applications (2nd ed.). New York: Routledge. 

Van Peet, A. A. J. (1992). De potentieeltheorie van intelligentie. [The potentiality theory of intelligence]. Amsterdam: University of Amsterdam, Ph.D. Thesis. 
}


\examples{
data("FamilyIQ")
str(FamilyIQ)
}
\keyword{datasets}
