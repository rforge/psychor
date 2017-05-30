\name{KoreanSpeech}
\alias{KoreanSpeech}
\docType{data}
\title{
Korean Speech Data
}
\description{
This dataset represents a subset of the data collected in an experiment on the phonetic profile of Korean formality by Winter and Grawunder (2012). The authors were interested in pitch changes between two different attitudes (formal vs. informal). 
}

\usage{
data("KoreanSpeech")

}

\format{
A data frame with 6 individuals (14 measurements per person) and the following variables: 
  \describe{
	\item{\code{subject}}{Subject ID}
	\item{\code{gender}}{Gender}
	\item{\code{scenario}}{7 interaction types ("making an appointment", "asking for a favor", "apologizing  for  coming  too  late", etc.)}
	\item{\code{attitude}}{Formality: formal vs. informal.}
	\item{\code{frequency}}{Pitch frequency in Hz}
}
}

\source{
Winter, B. (2013). Linear models and linear mixed effects models in R with linguistic applications. arXiv:1308.5499. (\url{http://arxiv.org/pdf/1308.5499.pdf}

Winter, B., & Grawunder, S. (2012) The phonetic profile of Korean formality. Journal of Phonetics, 40, 808-815. 
}


\examples{
data("KoreanSpeech")
str(KoreanSpeech)
}
\keyword{datasets}
