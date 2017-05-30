\name{Rmotivation2}
\alias{Rmotivation2}
\docType{data}
\title{
Psychometric structure of R package authors
}
\description{
This dataset contains factor scores (person parameters) based on a 2-PL IRT model fitted on the following three scales: word design questionnaire (WDQ; task, social, and knowledge characteristics), Reinholt's motivation scale (extrinsic, intrinsic, hyrbrid), and Schwartz' value scale (universalism, power, self-direction). 
}

\usage{
data("Rmotivation2")
}

\format{
A data frame with 764 individuals and the following 18 variables.
  \describe{
	\item{\code{lists}}{Participation in R lists.}
	\item{\code{meet}}{Participation in R meetings/conferences.}
	\item{\code{npkgs}}{Number of packages developed/contributed.}
	\item{\code{wtask}}{WDQ task subscale.}
	\item{\code{wsocial}}{WDQ social subscale.}
	\item{\code{wknowledge}}{WDQ knowledge subscale.}
	\item{\code{mextrinsic}}{Extrinsic motivation.}
	\item{\code{mhybrid}}{Hybrid motivation.}
	\item{\code{mintrinsic}}{Intrinsic motivation.}
	\item{\code{vuniversalism}}{Schwartz value universalism.}
	\item{\code{vpower}}{Schwartz value power.}
	\item{\code{vselfdirection}}{Schwartz value self-direction.}
	\item{\code{gender}}{Gender.}
	\item{\code{phd}}{PhD degree.}
	\item{\code{statseduc}}{Education in statistics.}
	\item{\code{fulltime}}{Full-time vs. part-time employment.}
	\item{\code{academia}}{Work in acedemia.}
	\item{\code{statswork}}{Work in the area of statistics.}
}
}

\source{
Mair, P., Hofmann, E., Gruber, K., Zeileis, A., & Hornik, K. (2015) Motivation, values, and work design as drivers of participation in the R open source Project for Statistical Computing. Proceedings of the National Academy of Sciences of the United States of America, 112(48), 14788-14792.
}

\seealso{\code{\link{Rmotivation}}}

\examples{
data("Rmotivation2")
str(Rmotivation2)
}
\keyword{datasets}
