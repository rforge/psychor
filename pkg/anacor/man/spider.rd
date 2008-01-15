\name{spider}
\alias{spider}
\docType{data}
\title{Hunting spider data}
\description{Abundance of hunting spiders in a Dutch dune area.
}
\usage{data(glass)}
\format{A list of data frames containing the frequency table (28 observations) and the row covariates.

Table:
  \describe{
  \item{\code{Alopacce}}{Abundance of Alopecosa accentuata.}
  \item{\code{Alopcune}}{Abundance of Alopecosa cuneata.}
  \item{\code{Alopfabr}}{Abundance of Alopecosa fabrilis.}
  \item{\code{Arctlute}}{Abundance of Arctosa lutetiana.}
  \item{\code{Arctperi}}{Abundance of Arctosa perita.}
  \item{\code{Auloalbi}}{Abundance of Aulonia albimana.}
  \item{\code{Pardlugu}}{Abundance of Pardosa lugubris.}
  \item{\code{Pardmont}}{Abundance of Pardosa monticola.}
  \item{\code{Pardnigr}}{Abundance of Pardosa nigriceps.}
  \item{\code{Pardpull}}{Abundance of Pardosa pullata.}
  \item{\code{Trocterr}}{Abundance of Trochosa terricola.}
  \item{\code{Zoraspin}}{Abundance of Zora spinimana.}
}

Row covariates:
 \describe{
  \item{\code{WaterCon}}{Log percentage of soil dry mass.}
  \item{\code{BareSand}}{Log percentage cover of bare sand.}
  \item{\code{FallTwig}}{Log percentage cover of fallen leaves and twigs.}
  \item{\code{CoveMoss}}{Log percentage cover of the moss layer.}
  \item{\code{CoveHerb}}{Log percentage cover of the herb layer.}
  \item{\code{ReflLux}}{Reflection of the soil surface with cloudless sky.}
  }
}

\references{Van der Aart, P.J.M. and Smeek-Enserink, N. (1975).
Correlations between distributions of hunting spiders (Lycosidae, Ctenidae)
and environmental characteristics in a dune area. Netherlands Journal of Zoology, 25, 1--45.
}
\examples{
data(spider)
str(spider)
}
\keyword{datasets}