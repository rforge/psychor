\name{tension}
\alias{tension}
\docType{data}
\title{
Perceived Tension in Music Over Time
}
\description{
This dataset comes from an experiment described Vines et al. (2006; the data were slightly modified). The authors were interested in how physical gestures of professional musicians contribute to the perception of emotion in a musical performance. 29 participants were exposed to the performance by either just listening (condition "auditory"), just seeing (condition "visual""), or both (condition "auditory-visual"). During the performance the participants had to move a slider to indicate the experienced tension they felt. They listened to the piece for 80 sec; every 10 msec the tension score (0 to 127) was recorded. This results in 800 tension measurement points per person (here provided as z-scores). 
}
\usage{
data("tension")
}

\format{
A data frame with 29 individuals and 800 measurement points. The last column \code{condition} contains the experimental conditions (auditory, visual, auditory-visual). 
}

\source{
Vines, B. W., Krumhansl, C. L., Wanderley, M. M., Levitin, D. J. (2006). Cross-modal interactions in the perception of musical performance. Cognition, 101, 80-113. 

Levitin, D. J., Nuzzo, R. L., Wines, B. W., & Ramsay, J. O. (2007). Introduction to functional data analysis. Canadian Psychology, 48, 135-155.
}

\examples{
data("tension")
str(tension)
}
\keyword{datasets}
