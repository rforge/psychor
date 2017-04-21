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
	\item{\code{Alcohol.1}}{Days with at least one drink in past year (T1).}
	\item{\code{Cigarettes.1}}{Number of cigarettes per day in past year (T1).}
	\item{\code{Alcohol2.1}}{Days with at least 5 drinks within a few hours in the past year (T1).}
	\item{\code{Marijuana.1}}{Times consumed marijuana in the past year (T1).}
	\item{\code{Alcohol.2}}{Days with at least one drink in past year (T2).}
	\item{\code{Cigarettes.2}}{Number of cigarettes per day in past year (T2).}
	\item{\code{Alcohol2.2}}{Days with at least 5 drinks within a few hours in the past year (T2).}
	\item{\code{Marijuana.2}}{Times consumed marijuana in the past year (T2).}
	\item{\code{Alcohol.3}}{Days with at least one drink in past year (T3).}
	\item{\code{Cigarettes.3}}{Number of cigarettes per day in past year (T3).}
	\item{\code{Alcohol2.3}}{Days with at least 5 drinks within a few hours in the past year (T3).}
	\item{\code{Marijuana.3}}{Times consumed marijuana in the past year (T3).}
	\item{\code{Alcohol.4}}{Days with at least one drink in past year (T4).}
	\item{\code{Cigarettes.4}}{Number of cigarettes per day in past year (T4).}
	\item{\code{Alcohol2.4}}{Days with at least 5 drinks within a few hours in the past year (T4).}
	\item{\code{Marijuana.4}}{Times consumed marijuana in the past year (T4).}
	\item{\code{Alcohol.5}}{Days with at least one drink in past year (T5).}
	\item{\code{Cigarettes.5}}{Number of cigarettes per day in past year (T5).}
	\item{\code{Alcohol2.5}}{Days with at least 5 drinks within a few hours in the past year (T5).}
	\item{\code{Marijuana.5}}{Times consumed marijuana in the past year (T5).}
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
