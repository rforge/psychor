\name{HRB}
\alias{HRB}
\docType{data}
\title{
Health Risk Behavior
}
\description{
Dataset based on a questionnaire assessing health risk behaviors, including smoking, drinking, and marijuana consumption. The questionnaire was presented to teenagers at 5 points in time (from middle school to high school). The items are binary: 0 = never, 
1 = at least one. 
}
\usage{
data("HRB")
}

\format{
A data frame with 538 individuals with 4 items presented at 5 points in time. Items: 
  \describe{
	\item{\code{Alcohol.1}}{Days with at least one drink in past year.}
	\item{\code{Cigarettes.1}}{Number of cigarettes per day in past year.}
	\item{\code{Alcohol2.1}}{Days with at least 5 drinks within a few hours in the past year.}
	\item{\code{Marijuana.1}}{Times consumed marijuana in the past year.}
}
}

\note{
Thanks to Peter Franz for providing this dataset.
}


\examples{
data("HRB")
str(HRB)
}
\keyword{datasets}
