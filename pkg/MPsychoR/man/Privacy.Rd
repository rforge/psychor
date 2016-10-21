\name{Privacy}
\alias{Privacy}
\docType{data}
\title{
Internet Privacy
}
\description{
These items measure various advantages and disadvantages which online users perceive when providing personal information on the Internet. The items are based on 25 qualitative interviews with online Marketing companies and experts as well as customer advocates. They represent the opinions of both organizations and individuals. Advantages of providing personal information online include support for purchasing decisions, increased satisfaction, targeted communication, participation in raffles, time savings and interesting content. Disadvantages include unsolicited 
advertising, excessive data collection, lack of information about data usage and decreasing service quality. 
}

\usage{
data("Privacy")
}

\format{
 A data frame with 405 individuals and the following 10 variables.
  \describe{
    \item{\code{apc1}}{Individualized communication supports me in making purchase decisions.}
    \item{\code{apc2}}{Individualized communication increases my satisfaction with the organization.}
    \item{\code{apc3}}{Individualization reduces the total amount of communication (e.g. the amount of emails I receive), since companies can advertise more target-oriented.}
    \item{\code{apc4}}{I provide correct data, if I have a change of winning prizes.}
    \item{\code{apc5}}{I provide correct data, if it saves me time (e.g. if I don't have to key in the data in the future).}
    \item{\code{apc6}}{I provide correct data, if I get access to interesting content.}
    \item{\code{dpc1}}{Participant's gender.}
    \item{\code{dpc2}}{On the Internet my data are permanently collected and I can do nothing against it.}
    \item{\code{dpc3}}{I feel that I am badly informed about the usage of my data.}
    \item{\code{dpc4}}{Increasing automation leads to a decrease in the service quality of companies.}
}
}

\source{
Treiblmaier, H. (2006) Datenqualitaet und individualisierte Kommunikation" [Data Quality and Individualized Communication], DUV Gabler Edition Wissenschaft, Wiesbaden.

Treiblmaier, H., Bentler, P. M., & Mair, P. (2011). Formative constructs implemented via common factors. Structural Equation Modeling: A Multidisciplinary Journal, 18, 1-17.
}


\examples{
data(Privacy)
str(Privacy)
}
\keyword{datasets}
