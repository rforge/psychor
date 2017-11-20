\name{Pashkam}
\alias{Pashkam}
\docType{data}
\title{
Goal-Directed Visual Processing
}
\description{
Data derived from an fMRI experiment on visual representations. In the original experiment there were three experimental conditions (color on objects and background, color on dots, color on objects), three brain regions of interest (V1, PFS, Superior IPS), and two tasks (color and shape). The data included here are two dissimilarity matrices involving eight objects presented to the participants. The first matrix is based on a color task, the second matrix on a shape task. 
}

\usage{
data("Pashkam")
}
\format{
  A list of 2 dissimilarity matrices (color task and shape task): 
  \describe{
	\item{\code{BD}}{Body}
	\item{\code{CT}}{Cat}
	\item{\code{CH}}{Chair}
	\item{\code{CR}}{Car}
	\item{\code{EL}}{Elephant}
	\item{\code{FA}}{Face}
	\item{\code{HO}}{House}
	\item{\code{SC}}{Scissors}
}}

\source{
Vaziri-Pashkam M., & Xu, Y. (2017) Goal-directed visual processing differentially impacts human ventral and dorsal visual representations. The Journal of Neuroscience, 37, 8767-8782.
}

\examples{
data(Pashkam)
str(Pashkam)
}
\keyword{datasets}
