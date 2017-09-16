\name{iatfaces}
\alias{iatfaces}
\docType{data}
\title{
Implicit Association Test (Faces)
}
\description{
The implicit association test (IAT) measures differential association of two target concepts with an attribute. In this experiment the participants saw images of people with long faces, images of people with wide faces, positively valenced words, and negatively valenced words. In the first critical block ("congruent block"), participants were asked to press one response key if they saw a long-faced person or a positive word and a different response key if they saw a wide-faced person or a negative word. In the second critical block ("incongruent block"), the pairing was reversed. Participants were asked to press one key for long-faced people or negative words and a different key for wide-faced people or positive words. IAT theory states that participants are expected to be able to respond fast in congruent conditions and slowly in incongruent conditions. The dataset contains trajectories of 4 participants. Each participant was exposed 80 trials: first, 40 congruent block trials, followed by 40 incongruent block trials. The response variable is latency. 
}

\usage{
data("iatfaces")
}

\format{
A data frame (4 individuals, 320 observations in total) with the following variables:
\describe{
    \item{\code{block}}{Congruent vs. incongruent.}
    \item{\code{latency}}{Response time latency.}
    \item{\code{id}}{Subject id.}
    \item{\code{trial}}{Trial number.}
}
}

\source{
Greenwald, A. G., & Banaji, M. R. (1995). Implicit social cognition: Attitudes, self-esteem, and stereotypes. Psychological Review, 102, 4-27.

Greenwald, A. G., McGhee, D.E., & Schwartz, J. K. L. (1998). Measuring individual differences in implicit cognition: The Implicit Association Test. Journal of Personality and Social Psychology, 74, 1464-1480.

Nosek, B. A., Banaji, M. R., & Greenwald, A. G. (2002). Harvesting implicit group attitudes and beliefs from a demonstration web site. Group Dynamics: Theory, Research, and Practice, 6, 101-115.
}

\note{
Thanks to Benedek Kurdi and Mahzarin Banaji for sharing this dataset. 
}

\examples{
data("iatfaces")
str(iatfaces)
}
\keyword{datasets}
