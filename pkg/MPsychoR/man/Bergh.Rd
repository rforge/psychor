\name{Bergh}
\alias{Bergh}
\docType{data}
\title{
Generalized Prejudice Dataset
}
\description{
Dataset from Bergh et al. (2016) where ethnic prejudice, sexism, sexual prejudice against gays and lesbians, and prejudice toward mentally people with disabilities are modeled as indicators of a generalized prejudice factor. It also includes indicators for agreeableness and openness. All variables are composite scores based on underlying 5-point questionnaire items. 
}

\usage{
data("Bergh")
}

\format{
A data frame with 861 individuals, 10 composite scores, and gender:
  \describe{
	\item{\code{EP}}{Ethnic prejudice}
	\item{\code{SP}}{Sexism}
	\item{\code{HP}}{Sexual prejudice against gays and lesbians}
	\item{\code{DP}}{Prejudice toward mentally people with disabilities}
	\item{\code{A1}}{Agreeableness indicator 1}
	\item{\code{A2}}{Agreeableness indicator 2}
	\item{\code{A3}}{Agreeableness indicator 3}
	\item{\code{O1}}{Openness indicator 1}
	\item{\code{O2}}{Openness indicator 2}
	\item{\code{O3}}{Openness indicator 3}
	\item{\code{gender}}{gender}
}
}

\source{
Bergh, R., Akrami, N., Sidanius, J., & Sibley, C. (2016) Is group membership necessary for understanding prejudice? A re-evaluation of generalized prejudice and its personality correlates. Journal of Personality and Social Psychology, 111, 367-395. 
}


\examples{
data("Bergh")
str(Bergh)
}
\keyword{datasets}
