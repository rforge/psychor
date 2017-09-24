\name{zareki}
\alias{zareki}
\docType{data}
\title{
Neuropsychological Test Battery for Number Processing and Calculation in Children
}
\description{
ZAREKI-R test battery (von Aster et al., 2006) for the assessment of dyscalculia in children. Includes subsets of 8 summation and 8 subtraction items, dichotomously scored, and 2 covariates.
}
\usage{data("zareki")}
\format{
  A data frame with 341 and 18 variables. Variables starting with \code{addit} are summation items, variables starting with \code{subtr} are subtraction items. \code{class} denotes elementary school class, \code{time} the time in min require to complete the test.
}

\source{
Koller, I., & Alexandrowicz, R. W. (2010) Eine psychometrische Analyse der ZAREKI-R mittels Rasch-Modellen [A psychometric analysis of the ZAREKI-R using Rasch-models]. Diagnostica 56, 57-67.

von Aster, M., Weinhold Zulauf, M., & Horn, R. (2006) Neuropsychologische Testbatterie fuer Zahlenverarbeitung und Rechnen bei Kindern (ZAREKI-R) [Neuropsychological Test Battery for Number Processing and Calculation in Children]. Harcourt Test Services, Frankfurt, Germany. 
}
\examples{
data(zareki)
str(zareki)
}
\keyword{datasets}
