\name{Rogers}
\alias{Rogers}
\docType{data}
\title{
Co-Morbid Obsessive-Compulsive Disorder and Depression
}
\description{
Depression/OCD Data Collected at Rogers Memorial Hospital. The scales used in this study were the Quick Inventory of Depressive Symptomatology - self-report version (QIDS-SR), and the Yale-Brown Obsessive Compulsive Scale - self-report (Y-BOCS-SR). The depression scale has 16 items (5 response categories), the OCD scale 10 items (4 response categories). 
}

\usage{
data("Rogers")
}

\format{
A data frame with 408 individuals and the following 26 variables (16 depression items followed by 10 OCD items).
  \describe{
	\item{\code{onset}}{Sleep-onset insomnia.}
	\item{\code{middle}}{Middle insomnia.}
	\item{\code{late}}{Early morning awakening.}
	\item{\code{hypersom}}{Hypersomnia.}
	\item{\code{sad}}{Sadness.}
	\item{\code{decappetite}}{Decreased appetite.}
	\item{\code{incappetite}}{Increased appetite.}
	\item{\code{weightloss}}{Weight loss.}
	\item{\code{weightgain}}{Weight gain.}
	\item{\code{concen}}{Concentration impairment.}
	\item{\code{guilt}}{Guilt and self-blame.}
	\item{\code{suicide}}{Suicidal thoughts, plans or attempts.}
	\item{\code{anhedonia}}{Anhedonia.}
	\item{\code{fatigue}}{Fatigue.}
	\item{\code{retard}}{Psychomotor retardation.}
	\item{\code{agitation}}{Agitation.}
	\item{\code{obtime}}{Time consumed by obsessions.}
	\item{\code{obinterfer}}{Interference due to obsessions.}
	\item{\code{obdistress}}{Distress caused by obsessions.}
	\item{\code{obresist}}{Difficulty resisting obsessions.}
	\item{\code{obcontrol}}{Difficulty controlling obsessions.}
	\item{\code{comptime}}{Time consumed by compulsions.}
	\item{\code{compinterf}}{Interference due to compulsions.}
	\item{\code{compdis}}{Distress caused by compulsions.}
	\item{\code{compresis}}{Difficulty resisting compulsions.}
	\item{\code{compcont}}{Difficulty controlling compulsions.}
	}
}

\source{
McNally, R. J., Mair, P., Mugno, B. L., and Riemann, B. C. (2017). Comorbid obsessive-compulsive disorder and depression: A Bayesian network approach. Psychological Medicine, 47(7), 1204-1214.
}


\examples{
data("Rogers")
str(Rogers)
}
\keyword{datasets}
