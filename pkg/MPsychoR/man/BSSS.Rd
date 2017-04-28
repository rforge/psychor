\name{BSSS}
\alias{BSSS}
\docType{data}
\title{
Brief Sensation Seeking Scale Questions (BSSS-8)
}
\description{
Haegeli et al. (2012) where interested in studying risk-taking behaviors of out-of-bounds skiers. The skiers where exposed to the ``Brief Sensation Seeking Scale'' (BSSS-8; Hoyle et al., 2002). It is a short 8-item scale with 5-point response categories. The scale has 4 subscales (with 2 items each): experience seeking (ES), boredom susceptibility (BS), thrill and adventure seeking (TAS) and disinhibition (DIS).
}

\usage{
data("BSSS")
}

\format{
A data frame with 1626 skiers and the following 8 items (5 response categories):
  \describe{
	 \item{\code{Explore}}{I would like to explore strange places.}
   \item{\code{Restless}}{I get restless when I spend too much time at home .}
	 \item{\code{Frightning}}{I like to do frightening things.}
	 \item{\code{Party}}{I like wild parties.}
	 \item{\code{Trip}}{I would like to take off on a trip with no pre-planned routes or timetables.}
	 \item{\code{Friends}}{I prefer friends who are exciting and unpredictable.}
	 \item{\code{Bungee}}{I would like to do bungee jumping.}
	 \item{\code{Illegal}}{I would love to have new and exciting experiences, even if they are illegal.}
}
}

\source{
Hoyle, R. H., Stephenson, M. T., Palmgreen, P., Lorch, E. P., & Donohew, R. L. (2002). Reliability and validity of a brief
measure of sensation. Personality and Individual Differences, 32, 401-414.

Haegeli, P., Gunn, M., & Haider, W. (2012). Identifying a high-risk cohort in a complex and dynamic risk environment: Out-of-bounds skiing--an example from avalanche safety. Prevention Science, 13, 562-573.
}


\examples{
data("BSSS")
str(BSSS)
}
\keyword{datasets}
