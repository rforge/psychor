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
     \item{\code{T0cdi1}}{I am sad all the time}
     \item{\code{T0cdi2r}}{Nothing will ever work out for me}
     \item{\code{T0cdi3}}{I do everything wrong}
     \item{\code{T0cdi4}}{Nothing is fun at all}
     \item{\code{T0cdi5r}}{I am bad all the time}
     \item{\code{T0cdi6}}{I am sure that terrible things will happen to me}
     \item{\code{T0cdi7r}}{I hate myself}
     \item{\code{T0cdi8r}}{All bad things are my fault}
     \item{\code{T0cdi10r}}{I feel like crying every day}
     \item{\code{T0cdi11r}}{Things bother me all the time}
     \item{\code{T0cdi12}}{I do not want to be with people at all}
     \item{\code{T0cdi13r}}{I cannot make up my mind about things}
     \item{\code{T0cdi14}}{I look ugly}
     \item{\code{T0cdi15r}}{I have to push myself all the time to do my schoolwork}
     \item{\code{T0cdi16r}}{I have trouble sleeping every night}
     \item{\code{T0cdi17}}{I am tired all the time}
     \item{\code{T0cdi18r}}{Most days I do not feel like eating}
     \item{\code{T0cdi19}}{I do not worry about aches and pains}
     \item{\code{T0cdi20}}{I do not feel alone}
     \item{\code{T0cdi21r}}{I never have fun at school}
     \item{\code{T0cdi22}}{I do not have any friends}
     \item{\code{T0cdi23}}{I do very badly in subjects I used to be good in}
     \item{\code{T0cdi24r}}{I can never be as good as other kids}
     \item{\code{T0cdi25r}}{Nobody really loves me}
     \item{\code{T0cdi26}}{I never do what I am told}
     \item{\code{T0cdi27}}{I get into fights all the time}
     \item{\code{race}}{Children's race} 
  }
}

\source{
Vaughn-Coaxum, R. A., Mair, P., & Weisz, J. R. (2015). Racial/ethnic differences in youth depression indicators: An Item Response Theory analysis of symptoms reported by White, Black, Asian, and Latino youths. Clinical Psychological Science, forthcoming. 
}
\examples{
data(YouthDep)
head(YouthDep)
str(YouthDep)
}
\keyword{datasets}
