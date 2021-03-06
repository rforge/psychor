\name{mcfa}
\alias{mcfa}
\title{ Two or more-sample CFA}
\description{
  Performs an analysis of configuration frequencies for two or more sets of counts. 
  \emph{This function is not designed to be called directly by the user but will only be used internally by \code{cfa()}.}
  Both the simple an the multiple cfa are handled by \code{cfa()}
}
\usage{
mcfa(cfg, cnts, sorton="chisq", sort.descending=TRUE, format.labels=TRUE)
}
\arguments{
  \item{cfg}{Contains the configurations. This can be a dataframe or a matrix. The dataframe can contain numbers,
              characters, factors or booleans. The matrix can consist of numbers, characters or booleans 
              (factors are implicitely re-converted to numerical levels). There must be >=3 columns.}
  \item{cnts}{ Contains the counts for the configuration. \code{cnts} is a matrix or dataframe with 2 or more columns.}
  \item{sorton}{ Determines the sorting order of the output. Can be set to \code{chisq}, \code{n}, or \code{label}.}
  \item{sort.descending}{Sort in descending order }
  \item{format.labels}{Format the labels of the configuration. This makes to output wider but it will increase the
        readability.}
}

\details{ This function is the "engine" \code{cfa()} will use. It does
the aggregation, summing up, and will calculate chi squared. All tests of significance are left to  \code{cfa()} 
}

\value{ The function returns the following list:
  \item{labels}{Configuration label} 
  \item{sums}{Sums for each configuration and each variable in the configuration}
  \item{counts}{Matrix of observed n of the given configuration}
  \item{expected}{Matrix of expected n for the given configuration}
  \item{chisq}{Chi squared for each configuration} 
}
\references{Krauth J., Lienert G. A. (1973, Reprint 1995) Die
     Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in
     Psychologie und Medizin, Beltz Psychologie Verlagsunion
     
     Lautsch, E., von Weber S. (1995) Methoden und Anwendungen der Konfigurationsfrequenzanalyse 
     in Psychologie und Medizin, Beltz Psychologie Verlagsunion

     Eye, A. von (1990) Introduction to configural frequency analysis.
     The search for types and anti-types in cross-classification. 
     Cambride 1990 }
     
\author{Stefan Funke <s.funke@t-online.de>}
\note{There are no hard-coded limits in the program so even large tables can be processed. 
     }

\seealso{ \code{\link{cfa}}, \code{\link{scfa}}}

\examples{ 
# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts1<-trunc(runif(250)*10) 
counts2<-trunc(runif(250)*10)
cfa(configs,cbind(counts1,counts2))
# cfa rather than mcfa!
}

\keyword{multivariate}
\keyword{htest}











