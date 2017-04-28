\name{granularity}
\alias{granularity}
\docType{data}
\title{
Granularity
}
\description{
Granularity refers to a person's ability to separate their emotions into specific types. People with low granularity struggle to separate their emotions (e.g., reporting that sadness, anger, fear, and others all just feel "bad""), whereas people with high granularity are very specific in how they parse their emotions (e.g., easily distinguishing between nuanced emotions like disappointment and frustration). A few outliers were removed compared to the original data.
}
\usage{
data("granularity")
}

\format{
 A data frame with 143 individuals and the following 3 variables. 
  \describe{
    \item{\code{gran}}{Granularity score}
    \item{\code{age}}{Participant's age}
    \item{\code{gender}}{Gender}
}
}


\examples{
data("granularity")
str(granularity)
}
\keyword{datasets}
