\name{ASTI}
\alias{ASTI}
\docType{data}
\title{
Adult Self-Transcendence Inventory
}
\description{
The ASTI (Levenson et al., 2005) is a self-report scale measuring the complex target construct of wisdom. The items can be assigned to five dimensions: self-knowledge and integration (SI), peace of mind (PM), non-attachment (NA), self-transcendence (ST), and presence in the here-and-now and growth (PG). 
}

\usage{
data("ASTI")
}

\format{
A data frame with 1215 individuals, 25 ASTI items (3 or 4 categories per items), and 2 covariates (gender, group). Item wordings: 
  \describe{
	\item{\code{ASTI1}}{I often engage in quiet contemplation. (PM; reversed)}
	\item{\code{ASTI2}}{I feel that my individual life is a part of a greater whole. (ST)}
	\item{\code{ASTI3}}{I don't worry about other people's opinions of me. (NA)}
	\item{\code{ASTI4}}{I feel a sense of belonging with both earlier and future generations. (ST)}
	\item{\code{ASTI5}}{My peace of mind is not easily upset. (PM)}
	\item{\code{ASTI6}}{My sense of well-being does not depend on a busy social life. (NA)}
	\item{\code{ASTI7}}{I feel part of something greater than myself. (ST)}
	\item{\code{ASTI8}}{My happiness is not dependent on other people and things. (NA; reversed)}
	\item{\code{ASTI9}}{I do not become angry easily. (PM)}
	\item{\code{ASTI10}}{I have a good sense of humor about myself. (SI; reversed)}
	\item{\code{ASTI11}}{I find much joy in life. (PG; reversed)}
	\item{\code{ASTI12}}{Material possessions don't mean much to me. (NA)}
	\item{\code{ASTI13}}{I feel compassionate even toward people who have been unkind to me. (ST)}
	\item{\code{ASTI14}}{I am not often fearful. (PG)}
	\item{\code{ASTI15}}{I can learn a lot from others. (PG)}
	\item{\code{ASTI16}}{I often have a sense of oneness with nature. (ST)}
	\item{\code{ASTI17}}{I am able to accept my mortality. (PG)}
	\item{\code{ASTI18}}{I often "lose myself" in what I am doing. (PG)}
	\item{\code{ASTI19}}{I feel that I know myself. (SI; reversed)}
	\item{\code{ASTI20}}{I am accepting of myself, including my faults. (SI; reversed)}
	\item{\code{ASTI21}}{I am able to integrate the different aspects of my life. (SI; reversed)}
	\item{\code{ASTI22}}{I can accept the impermanence of things. (PM; reversed)}
	\item{\code{ASTI23}}{I have grown as a result of losses I have suffered. (PG; reversed)}
	\item{\code{ASTI24}}{Whatever [good] I do for others, I do for myself. (ST; reversed)}
	\item{\code{ASTI25}}{Whatever [bad] I do to others, I do to myself. (ST)}
	\item{\code{gender}}{gender}
	\item{\code{group}}{student vs. non-student}
}
}

\source{
Levenson, M. R., Jennings, P. A., Aldwin, C. M., & Shiraishi, R. W. (2005). Self-transcendence: conceptualization and measurement. The International Journal of Aging and Human Development, 60, 127-143. 

Koller I., Levenson, M. R. , & Glueck, J. (2017). What do you think you are measuring? A mixed-methods procedure for assessing the content validity of test items and theory-based scaling. Frontiers in Psychology, 8(126), 1-20. 
}


\examples{
data(ASTI)
si <- ASTI[ ,c(10,19,20,21)]            ## self-knowledge and integration
pm <- ASTI[ ,c(1,5,9,22)]               ## peace of mind
na <- ASTI[ ,c(3,6,8,12)]               ## non-attachment
st <- ASTI[ ,c(2,4,7,13,16,24,25)]      ## self-transcendence
pg <- ASTI[ ,c(11,14,15,17,18,23)]      ## Presence in the here-and-now and growth
}
\keyword{datasets}
