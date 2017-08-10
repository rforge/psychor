\name{Lakes}
\alias{Lakes}
\docType{data}
\title{
Response to challenge scale
}
\description{
The response to challenge scale (RCS) is a theory-derived, observer-rated measure of children's self-regulation in response to a physically challenging situation (Lakes & Hoyt, 2004; Lakes, 2012). It asks raters to make inferences in 3 domains: cognitive (6 items), affective/motivational (7 items), and physical (3 items). The data included here are post test ratings from the study presented in Lakes & Hoyt (2009). 
}

\usage{
data("Lakes")
}

\format{
A data frame in long format with 194 individuals and the following variables:
  \describe{
	\item{\code{personID}}{Person ID.}
	\item{\code{raterID}}{Rater ID.}
	\item{\code{item}}{Items for 3 subtests.}
	\item{\code{score}}{7-point response score.}
	\item{\code{subtest}}{Subtests (cognitive, affective, physical).}
}
}

\source{
Lakes, K. D. (2012). The Response to Challenge Scale (RCS): The development and construct validity of an observer-rated measure of children's 
self-regulation. The International Journal of Educational and Psychological Assessment, 10, 83-96. 

Lakes, K. D, & Hoyt, W. T. (2004). Promoting self-regulation through school-based martial arts training. Journal of Applied Developmental Psychology, 25, 283-302. 

Lakes, K. D., & Hoyt, W. T. (2009). Applications of generalizability theory to clinical child and adolescent psychology research. Journal of Clinical Child & Adolescent Psychology, 38, 144-165. 
}


\examples{
data("Lakes")
str(Lakes)
}
\keyword{datasets}
