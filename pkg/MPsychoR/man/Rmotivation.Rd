\name{Rmotivation}
\alias{Rmotivation}
\docType{data}
\title{
Motivational structure of R package authors
}
\description{
Motivation is accurately understood as a complex continuum of intrinsic, extrinsic, and internalized extrinsic motives. This dataset 
contains three subscales for that measure extrinsic (12 items), hybrid (19 items), and intrinsic (5 items) aspects of motivation 
in relation to why package authors contribute to the R environment. The items were taken from Reinholt's motivation scale and adapted 
to R package authors. Each item started with "I develop R packages, because...".
}

\usage{
data("Rmotivation")
}

\format{
A data frame with 852 individuals, 36 motivation items, and 9 covariates:
  \describe{
	\item{\code{ext1}}{I can publish the packages in scientific journals.}
	\item{\code{ext2}}{they are part of my master / PhD thesis.}
	\item{\code{ext3}}{I need them for teaching courses.}
	\item{\code{ext4}}{I develop them for clients who pay me.}
	\item{\code{ext5}}{they are a byproduct of my empirical research. If I cannot find suitable existing software to analyze my data, I develop software components myself.}
	\item{\code{ext6}}{they are a byproduct of my methodological research. If I develop/extend methods, I develop accompanying software, e.g., for illustrations and simulations.}
	\item{\code{ext7}}{I expect an enhancement of my career from it.}
	\item{\code{ext8}}{my employer pays me to do so.}
	\item{\code{ext9}}{that's what my friends do.}
	\item{\code{ext10}}{it is expected from me.}
	\item{\code{ext11}}{that's what my work colleagues do.}
	\item{\code{ext12}}{it comes more or less with my job.}
	\item{\code{hyb1}}{it is an important task for me.}
	\item{\code{hyb2}}{I believe that it is a necessity.}
	\item{\code{hyb3}}{I believe it is vital to improve R.}
	\item{\code{hyb4}}{I feel that R requires continuous enhancement.}
	\item{\code{hyb5}}{I think that it is of importance.}
	\item{\code{hyb6}}{it is part of my identity.}
	\item{\code{hyb7}}{it is important for my personal goals but for no apparent rewards, such as money, career opportunities, etc.}
	\item{\code{hyb8}}{it is part of my character to do so.}
	\item{\code{hyb9}}{it is an integral part of my personality.}
	\item{\code{hyb10}}{it is in line with my personal values.}
	\item{\code{hyb11}}{I feel an obligation towards the R community.}
	\item{\code{hyb12}}{it reflects my responsibility towards the R community.}
	\item{\code{hyb13}}{I believe that it is appropriate to do so.}
	\item{\code{hyb14}}{I aim for social approval of my activities.}
	\item{\code{hyb15}}{I am committed to the R community.}
	\item{\code{hyb16}}{I can feel satisfied with my performance.}
	\item{\code{hyb17}}{it leaves me with a feeling of accomplishment.}
	\item{\code{hyb18}}{it gives me satisfaction to produce something of high quality.}
	\item{\code{hyb19}}{I get the feeling that I've accomplished something of great value.}
	\item{\code{int1}}{I enjoy undertaking the required tasks.}
	\item{\code{int2}}{I take pleasure in applying my skills.}
	\item{\code{int3}}{it means pure fun for me.}
	\item{\code{int4}}{I feel that it is an interesting exercise.}
	\item{\code{int5}}{it is a joyful activity.}
	\item{\code{lists}}{Participation in R lists.}
	\item{\code{meet}}{Participation in R meetings/conferences.}
	\item{\code{npkgs}}{Number of packages developed/contributed.}
	\item{\code{gender}}{Gender.}
	\item{\code{phd}}{PhD degree.}
	\item{\code{statseduc}}{Education in statistics.}
	\item{\code{fulltime}}{Full-time vs. part-time employment.}
	\item{\code{academia}}{Work in acedemia.}
	\item{\code{statswork}}{Work in the area of statistics.}
}
}

\source{
Mair, P., Hofmann, E., Gruber, K., Zeileis, A., & Hornik, K. (2015) Motivation, values, and work design as drivers of participation in the R open source Project for Statistical Computing. Proceedings of the National Academy of Sciences of the United States of America, 112(48), 14788-14792.

Reinholt, M. (2006). No more polarization, please! Towards a more nuanced perspective on motivation in organizations. Technical report, Center for Strategic Management Working Paper Series, Copenhagen Business School, Copenhagen, Denmark.
}


\examples{
data(Rmotivation)
str(Rmotivation)
}
\keyword{datasets}
