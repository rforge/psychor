\name{YouthDep}
\alias{YouthDep}
\docType{data}
\title{
Youth Depression Indicators
}
\description{
Contains Children's Depression Inventory (CDI) measures of sixth and seventh grade students. In total, there are 26 CDI items 
(on of the original CDI items asking about suicidal ideation was removed) with three response categories each
(e.g., 0 = nobody really loves me, 1 = I am not sure if anybody loves me, or 2 = I am sure that somebody loves me).
}
\usage{data("YouthDep")}
\format{
  A data frame with 2290 on the following 27 variables.
  \describe{
     \item{\code{CDI1}}{I am sad all the time}
     \item{\code{CDI2r}}{Nothing will ever work out for me}
     \item{\code{CDI3}}{I do everything wrong}
     \item{\code{CDI4}}{Nothing is fun at all}
     \item{\code{CDI5r}}{I am bad all the time}
     \item{\code{CDI6}}{I am sure that terrible things will happen to me}
     \item{\code{CDI7r}}{I hate myself}
     \item{\code{CDI8r}}{All bad things are my fault}
     \item{\code{CDI10r}}{I feel like crying every day}
     \item{\code{CDI11r}}{Things bother me all the time}
     \item{\code{CDI12}}{I do not want to be with people at all}
     \item{\code{CDI13r}}{I cannot make up my mind about things}
     \item{\code{CDI14}}{I look ugly}
     \item{\code{CDI15r}}{I have to push myself all the time to do my schoolwork}
     \item{\code{CDI16r}}{I have trouble sleeping every night}
     \item{\code{CDI17}}{I am tired all the time}
     \item{\code{CDI18r}}{Most days I do not feel like eating}
     \item{\code{CDI19}}{I do not worry about aches and pains}
     \item{\code{CDI20}}{I do not feel alone}
     \item{\code{CDI21r}}{I never have fun at school}
     \item{\code{CDI22}}{I do not have any friends}
     \item{\code{CDI23}}{I do very badly in subjects I used to be good in}
     \item{\code{CDI24r}}{I can never be as good as other kids}
     \item{\code{CDI25r}}{Nobody really loves me}
     \item{\code{CDI26}}{I never do what I am told}
     \item{\code{CDI27}}{I get into fights all the time}
     \item{\code{race}}{Children's race} 
  }
}

\source{
Vaughn-Coaxum, R. A., Mair, P., & Weisz, J. R. (2015). Racial/ethnic differences in youth depression indicators: An Item Response Theory analysis of symptoms reported by White, Black, Asian, and Latino youths. Clinical Psychological Science, 4, 239-253. 
}
\examples{
data(YouthDep)
head(YouthDep)
str(YouthDep)
}
\keyword{datasets}
