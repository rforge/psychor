\name{RWDQ}
\alias{RWDQ}
\docType{data}
\title{
Work design questionnaire R package authors
}
\description{
Contains the knowledge characteristics subscale of the Work Design Questionnaire (WDQ). Knowledge characteristics include job complexity, information processing, problem solving, skill variety, and specialization. 
}

\usage{
data("RWDQ")
}

\format{
A data frame with 1055 individuals and 18 items: job complexity (22-24), information processing (25-27), problem solving (28-31), variety of skills (32-35), specialization (36-39). Item wordings: 
  \describe{
	\item{\code{wdq_22}}{The work on R packages requires that I only do one task or activity at a time.}
	\item{\code{wdq_23}}{The work on R packages comprises relatively uncomplicated tasks.}
	\item{\code{wdq_24}}{The work on R packages involves performing relatively simple tasks.}
	\item{\code{wdq_25}}{The work on R packages requires that I engage in a large amount of thinking.}
	\item{\code{wdq_26}}{The work on R packages requires me to keep track of more than one thing at a time.}
	\item{\code{wdq_27}}{The work on R packages requires me to analyze a lot of information}
	\item{\code{wdq_28}}{The work on R packages involves solving problems that have no obvious correct answer.}
	\item{\code{wdq_29}}{The work on R packages requires me to be creative.}
	\item{\code{wdq_30}}{The work on R packages often involves dealing with problems that I have not encountered before.}
	\item{\code{wdq_31}}{The work on R packages requires unique ideas or solutions to problems.}
	\item{\code{wdq_32}}{The work on R packages requires data analysis skills.}
	\item{\code{wdq_33}}{The work on R packages requires programming skills.}
	\item{\code{wdq_34}}{The work on R packages requires technical skills regarding package building and documentation.}
	\item{\code{wdq_35}}{The work on R packages requires the use of a number of skills.}
	\item{\code{wdq_36}}{The work on R packages is highly specialized in terms of purpose, tasks, or activities.}
	\item{\code{wdq_37}}{The tools, procedures, materials, and so forth used to develop R packages are highly specialized in terms of purpose.}
	\item{\code{wdq_38}}{The work on R packages requires very specialized knowledge.}
	\item{\code{wdq_39}}{The work on R packages requires a depth of expertise.}
}
}

\source{
Mair, P., Hofmann, E., Gruber, K., Zeileis, A., & Hornik, K. (2015) Motivation, values, and work design as drivers of participation in the R open source Project for Statistical Computing. Proceedings of the National Academy of Sciences of the United States of America, 112(48), 14788-14792.

Morgeson, F. P., & Humphrey, S. E. (2006). The Work Design Questionnaire (WDQ): Developing and validating a comprehensive measure for assessing job design and the nature of work. Journal of Applied Psychology, 91, 1321-1339
}


\examples{
data(RWDQ)
str(RWDQ)
}
\keyword{datasets}
