\name{ageiat}
\alias{ageiat}
\docType{data}
\title{
Time Series Implicit Association Test (Age)
}
\description{
The implicit association test (IAT) measures differential association of two target concepts with an attribute. The outcome measure is the IAT score D (similar to Cohen's d effect size). There are different types of IAT. This dataset contains D-measures from the Age IAT (where most individuals have an implicit preference for young over old). They were collected on the ProjectImplicit platform (\url{http://implicit.harvard.edu/}) from May 2005 to December 2015. Within each each month the participants D-measures were averaged. This leads to a time series with 140 observations.
}

\usage{
data("ageiat")
}

\format{
A vector of D-scores, measured at 140 points in time (May 2005 - December 2015). 
}

\source{
Greenwald, A. G., & Banaji, M. R. (1995). Implicit social cognition: Attitudes, self-esteem, and stereotypes. Psychological Review, 102, 4-27.

Greenwald, A. G., McGhee, D.E., & Schwartz, J. K. L. (1998). Measuring individual differences in implicit cognition: The Implicit Association Test. Journal of Personality and Social Psychology, 74, 1464-1480.

Nosek, B. A., Banaji, M. R., & Greenwald, A. G. (2002). Harvesting implicit group attitudes and beliefs from a demonstration web site. Group Dynamics: Theory, Research, and Practice, 6, 101-115.
}

\note{
Thanks to Tessa Charlesworth and Mahzarin Banaji for sharing this dataset. 
}

\examples{
data("ageiat")
str(ageiat)
}
\keyword{datasets}
