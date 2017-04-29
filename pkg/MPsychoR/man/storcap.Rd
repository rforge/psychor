\name{storcap}
\alias{storcap}
\docType{data}
\title{
EEG Visual Working Memory Storage Capacity
}
\description{
The data were collected in an experiment on visual working memory storage capacity. The left-right electrode voltage were averaged. The sampling frequency was 2 Hz. There were 6 conditions in the experiment: Set Size 1 - Ipsilateral Activity; Set Size 1 - Contralateral Activity; Set Size 3 - Ipsilateral Activity; Set Size 3 - Contralateral Activity; Set Size 1 - Contralateral Delay Activity; Set Size 3 - Contralateral Delay Activity. Memory display from 0-300 msec, consolidation period 300-1200 msec, after 1200 msec test period.
}
\usage{
data("storcap")
}

\format{
 A data frame containing the following variables
  \describe{
    \item{\code{id}}{Subject ID}
    \item{\code{channel}}{EEG channel (13 in total)}
    \item{\code{time}}{Time}
    \item{\code{cond}}{Experimental conditions}
    \item{\code{voltage}}{Voltage electrode}
}
}

\note{
Thanks to Hrag Pailian for providing this dataset.
}


\examples{
data("storcap")
str(storcap)
}
\keyword{datasets}
