\name{learnemo}
\alias{learnemo}
\docType{data}
\title{
Learning related emotions in mathematics
}
\description{
This dataset considers achievement emotions students typically experience when learning mathematics. The authors considered 5 emotions: enjoyment (coded as 1), pride (2), anger (3), anxiety (4) and boredom (5). The data are organized in terms of paired comparisons (in standard order).
}

\usage{
data("learnemo")
}

\format{
  A data frame with 111 individuals and the following paired comparisons (0 if the first emotion was chosen, 2 if the second emotion was chosen, and 1 if no decision was made).
  \describe{
	\item{\code{pc1_2}}{enjoyment vs. pride.}
	\item{\code{pc1_3}}{enjoyment vs. anger.}
	\item{\code{pc2_3}}{pride vs. anger.}
	\item{\code{pc1_4}}{enjoyment vs. anxiety.}
	\item{\code{pc2_4}}{pride vs. anxiety.}
	\item{\code{pc3_4}}{anger vs. anxiety.}
	\item{\code{pc1_5}}{enjoyment vs. boredom.}
	\item{\code{pc2_5}}{pride vs. boredom.}
	\item{\code{pc3_5}}{anger vs. boredom.}
	\item{\code{pc4_5}}{anxiety vs. boredom.}
	\item{\code{sex}}{Participant's sex (1 = male, 2 = female).}
}
}

\source{
Grand, A., & Dittrich, R. (2015) Modelling assumed metric paired comparison data - application to learning related emotions. Austrian Journal of Statistics, 44, 3-15. 
}


\examples{
data("learnemo")
str(learnemo)
}
\keyword{datasets}
