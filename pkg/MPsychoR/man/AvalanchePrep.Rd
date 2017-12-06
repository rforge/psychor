\name{AvalanchePrep}
\alias{AvalanchePrep}
\docType{data}
\title{
Preparanedness Backcountry Skiing
}
\description{
Haegeli et al. (2012) studied high-risk cohorts in a complex and dynamic risk environment. This dataset contains four variables related to preparedness before going backcountry skiing. The variables with response categories are are 1) check avalanche danger information (check conditions on internet prior to leaving home; talk to ski patrol; check postings at gates or information kiosks at resort; do not check or Do not know), 2) discuss avalanche hazard in your group (all the time; 50\% to 90\% of time; 10\% to 40\% of time; never or solo traveller), 3) approach to decision making (dedicated leader or everybody contributes; person in front decides; everybody makes their own choices or solo traveller), and 4) use of avalanche safety gear (everybody carries beacon, shovel and probe; everybody carries beacon or beacon and shovel; some in group carry beacons; some in group have cell phones; no safety equipment is carried).
}

\usage{
data("AvalanchePrep")
}

\format{
A data frame with 1355 skiers and the following 4 items:
  \describe{
	 \item{\code{info}}{Check avalanche danger information.}
   \item{\code{discuss}}{Discuss avalanche hazard in your group.}
	 \item{\code{gear}}{Use of avalanche safety gear.}
	 \item{\code{decision}}{Approach to decision making.}
}
}

\source{
Haegeli, P., Gunn, M., & Haider, W. (2012). Identifying a high-risk cohort in a complex and dynamic risk environment: Out-of-bounds skiing--an example from avalanche safety. Prevention Science, 13, 562-573.
}


\examples{
data("AvalanchePrep")
str(AvalanchePrep)
}
\keyword{datasets}
