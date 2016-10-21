\name{BrainIQ}
\alias{BrainIQ}
\docType{data}
\title{
Brain Size and Intelligence
}
\description{
Willerman et al. (1991) conducted their study at a large southwestern university. They selected a sample of 40 right-handed Anglo introductory psychology students who had indicated no history of alcoholism, unconsciousness, brain damage, epilepsy, or heart disease. These subjects were drawn from a larger pool of introductory psychology students with total Scholastic Aptitude Test Scores higher than 1350 or lower than 940 who had agreed to satisfy a course requirement by allowing the administration of four subtests (Vocabulary, Similarities, Block
Design, and Picture Completion) of the Wechsler (1981) Adult Intelligence Scale-Revised. With prior approval of the University's research review board, students selected for MRI were required to obtain prorated full-scale IQs of greater than 130 or less than 103, and were equally divided by sex and IQ classification.

}
\usage{
data("BrainIQ")
}

\format{
 A data frame with 40 individuals and the following 7 variables.
  \describe{
    \item{\code{Gender}}{Participant's gender.}
    \item{\code{FSIQ}}{Full Scale IQ.}
    \item{\code{VIQ}}{Verbal IQ.}
    \item{\code{PIQ}}{Performance IQ.}
    \item{\code{Weight}}{Body weight.}
    \item{\code{Height}}{Body height.}
    \item{\code{MRI_Count}}{MRI pixel count (brain size).}
}
}

\source{
Willerman, L., Schultz, R., Rutledge, J. N., & Bigler, E. (1991). In vivo brain size and intelligence. Intelligence, 15, 223-228.
}


\examples{
data(BrainIQ)
str(BrainIQ)
}
\keyword{datasets}
