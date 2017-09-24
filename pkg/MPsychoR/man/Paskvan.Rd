\name{Paskvan}
\alias{Paskvan}
\docType{data}
\title{
Cognitive appraisal of work intensification
}
\description{
Due to economic and technological changes, work has intensified over the past few decades. This intensification of work takes a toll on 
employees well-being and job satisfaction. Paskvan et al. (2016) established a model which explores the effects of work intensification on various outcomes (emotional exhaustion, job satisfaction). They used cognitive appraisal (i.e., how an individual views a situation) as a mediator and the participative climate as a moderator of the relationship between work intensification and cognitive appraisal. 
}

\usage{
data("Paskvan")
}

\format{
A data frame with 803 individuals and the following 4 variables.
  \describe{
	\item{\code{pclimate}}{Participative climate.}
	\item{\code{wintense}}{Work Intensification.}
	\item{\code{cogapp}}{Cognitive appraisal of work intensification.}
	\item{\code{emotion}}{Emotional exhaustion.}
	}
}

\source{
Paskvan, M., Kubicek, B., Prem, R., & Korunka, C. (2016). Cognitive appraisal of work intensification. International Journal of Stress Management, 23, 124-146.
}


\examples{
data("Paskvan")
str(Paskvan)
}
\keyword{datasets}
