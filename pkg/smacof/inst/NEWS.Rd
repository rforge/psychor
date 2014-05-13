\name{NEWS}
\title{News for Package \pkg{smacof}}
\encoding{UTF-8}
\section{Changes in Version 1.4-0}{
  \itemize{
    \item Gravity model for dissimilarity computation based on co-occurrences added: see \code{gravity()}. 
    \item The \code{circle} argument in \code{smacofRect()} allows for circular restrictions on the configurations.
    \item All smacof functions now return Stress-1 (before it was squared Stress-1); the value return by the function is now called \code{stress}. 
    \item The \code{metric} argument is deprecated. Now \code{type} should be used (either \code{"ordinal"}, \code{"interval"}, or \code{"ratio"})
    for pretty much all smacof functions. 
    \item \code{smacofSphere.primal()} and \code{smacofSphere.dual()} are now merged into a single function \code{smacofSphere()} with a corresponding \code{algorithm} argument. 
    \item \code{smacofConstraint()} now returns weight matrix \code{C}.
    \item Plot updates:
    \itemize{
       \item For unfolding (\code{smacofRect()}) configurations: \code{joint = TRUE} as default. 
       \item Configuration plot looks a bit fancier now (labels and points as default), also aspect ratio = 1 as default.
       \item A note in plot help file added that in the stress plot the SPP from a various smacof functions is rescaled to percentages (which sum up to 100). 
    }
    \item Datasets:
    \itemize{
       \item Labels for rectangle dataset now provided in English.
       \item Dataset on Republican statements added (document-term matrix).  
    }
  }
}
