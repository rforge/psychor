\name{CEAQ}
\alias{CEAQ}
\docType{data}
\title{
Children's Empathic Attitudes Questionnaire (CEAQ)
}
\description{
The CEAQ (Funk et al., 2008) is a scale to measure empathy of late elementary and middle-school aged children.
}

\usage{
data("CEAQ")
}

\format{
A data frame with 208 children, 16 CEAQ items and 3 covariates (age, grade, gender): Item wordings: 
  \describe{
  \item{\code{ceaq1}}{When I'm mean to someone, I usually feel bad about it later.}
  \item{\code{ceaq2}}{I'm happy when the teacher says my friend did a good job.}
	\item{\code{ceaq3}}{I would get upset if I saw someone hurt an animal.}
	\item{\code{ceaq4}}{I understand how other kids feel.}
	\item{\code{ceaq5}}{I would feel bad if my mom's friend got sick.}
	\item{\code{ceaq6}}{Other people's problems really bother me.}
	\item{\code{ceaq7}}{I feel happy when my friend gets a good grade.}
	\item{\code{ceaq8}}{When I see a kid who is upset it really bothers me.}
	\item{\code{ceaq9}}{I would feel bad if the kid sitting next to me got in trouble.}
	\item{\code{ceaq10}}{It's easy for me to tell when my mom or dad has a good day at work.}
	\item{\code{ceaq11}}{It bothers me when my teacher doesn't feel well.}
	\item{\code{ceaq12}}{I feel sorry for kids who can't find anyone to hang out with.}
	\item{\code{ceaq13}}{Seeing a kid who is crying makes me feel like crying.}
	\item{\code{ceaq14}}{If two kids are fighting, someone should stop it.}
	\item{\code{ceaq15}}{It would bother me if my friend got grounded.}
	\item{\code{ceaq16}}{When I see someone who is happy, I feel happy too.}
	\item{\code{age}}{Children's age.}
	\item{\code{grade}}{Children's grade.}
	\item{\code{gender}}{Gender.}
}
}

\source{
Funk, J. B., Fox, C. M., Chang, M., & Curtiss, K. (2008). The development of the Children's Empathic Attitudes Questionnaire using classical and Rasch analyses. Journal of Applied Developmental Psychology, 29, 187-196. 

Bond, T. G., & Fox, C. M. (2015). Applying the Rasch Model: Fundamental Measurement in the Human Sciences. Routledge.
}


\examples{
data(CEAQ)
str(CEAQ)
}
\keyword{datasets}
