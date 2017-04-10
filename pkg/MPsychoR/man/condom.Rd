\name{condom}
\alias{condom}
\docType{data}
\title{
Attitude towards condoms
}
\description{
This dataset is a modified version of the data used in de Ayala (2009). Originally, the data come from the voluntary HIV counseling and testing efficacy study performed by the center for AIDS prevention studies (2003). 
}
\usage{
data("condom")
}

\format{
 A data frame with 500 individuals and the following 7 variables. The 6 items were scored on a 4-point response scale 
 (0 ... strongly disagree; 4 ... strongly agree).
  \describe{
    \item{\code{Feel}}{Condom does not have a good feel.}
    \item{\code{Buy}}{I am embarrassed to buy condoms.}
    \item{\code{Put}}{I am embarrased to put on condom.}
    \item{\code{Break}}{Condoms break/slip off.}
    \item{\code{Cheat}}{My partner wants condoms to cheat.}
    \item{\code{Uncomfortable}}{My friends said that condoms are uncomfortable.}
    \item{\code{Country}}{Participant's country (artificially added).}
}
}

\source{
de Ayala, R. J. (2009). The Theory and Practice of Item Response Theory. Guilford Press, New York
}


\examples{
data(condom)
str(condom)
}
\keyword{datasets}
