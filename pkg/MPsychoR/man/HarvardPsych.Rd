\name{HarvardPsych}
\alias{HarvardPsych}
\docType{data}
\title{
Research Topics Harvard Psychology Faculty
}
\description{
A frequency table with the faculty members in the rows and the research topics in the colunms. The data are based on a scraping job 
from the faculty website by extracting the research summary of each faculty members. Subsequently, the data were cleaned using 
basic text processing tools. Finally, a document term matrix was created containing the most important keywords in the columns.
}
\usage{
data("HarvardPsych")
}
\format{
 A word frequency table spanned 29 faculty members and 43 keywords.
}

\source{
URL: \url{http://psychology.fas.harvard.edu/faculty}
}

\examples{
data("HarvardPsych")
str(HarvardPsych)
}
\keyword{datasets}
