\name{Wenchuan}
\alias{Wenchuan}
\docType{data}
\title{
Wenchuan PTSD Dataset
}
\description{
PTSD (posttraumatic stress disorder) symptoms reported by survivors of the Wenchuan earthquake in China using the PTSD checklist-civilian version (PCL-C). All items were scaled on a 5-point Likert scale (1 ... not at all; 2 ... a little bit; 3 ... moderately; 4 ... quite a bit;
5 ... extremely).
}
\usage{data("Wenchuan")}
\format{
  A data frame with 362 observations on the following 17 variables.
  \describe{
    \item{\code{intrusion}}{Repeated, disturbing memories, thoughts, or images of a stressful experience from the past?}
    \item{\code{dreams}}{Repeated, disturbing dreams of a stressful experience from the past?}
    \item{\code{flash}}{Suddenly acting or feeling as if a stressful experience were happening again (as if you were reliving it)?}
    \item{\code{upset}}{Feeling very upset when something reminded you of a stressful experience from the past?}
    \item{\code{physior}}{Having physical reactions (e.g., heart pounding, trouble breathing, sweating) when something reminded you of a stressful experience from the past? }
    \item{\code{avoidth}}{Avoiding thinking about or talking about a stressful experience from the past or avoiding having feelings related to it?}
    \item{\code{avoidact}}{Avoiding activities or situations because they reminded you of a stressful experience from the past?}
    \item{\code{amnesia}}{Trouble remembering important parts of a stressful experience from the past?}
    \item{\code{lossint}}{Loss of interest in activities that you used to enjoy?}
    \item{\code{distant}}{Feeling distant or cut off from other people?}
    \item{\code{numb}}{Feeling emotionally numb or being unable to have loving feelings for those close to you?}
    \item{\code{future}}{Feeling as if your future will somehow be cut short?}
    \item{\code{sleep}}{Trouble falling or staying asleep?}
    \item{\code{anger}}{Feeling irritable or having angry outbursts?}
    \item{\code{concen}}{Having difficulty concentrating?}
    \item{\code{hyper}}{Being "super-alert" or watchful or on guard?}
    \item{\code{startle}}{Feeling jumpy or easily startled?}
  }
}

\source{
McNally, R. J., Robinaugh, D. J., Wu, G. W. Y., Wang, L., Deserno, M. K., & Borsboom, D. (2015). Mental disorders as causal systems: A network approach to posttraumatic stress disorder. Clinical Psychological Science, forthcoming. 
}
\examples{
data(Wenchuan)
head(Wenchuan)
str(Wenchuan)
}
\keyword{datasets}
