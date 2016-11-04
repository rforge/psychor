\name{NeuralActivity}
\alias{NeuralActivity}
\alias{NeuralScales}
\alias{NeuralScanner}
\alias{NeuralHM}
\docType{data}
\title{
Neural Activity
}
\description{
20 participants were scanned (fMRI) while performing a task designed to elicit their thoughts about 60 mental states. 
On each trial, participants saw the name of a mental state (e.g., "awe"), and decided which of two scenarios would better evoke that mental state in another person (e.g., "seeing the Pyramids" or "watching a meteor shower"). Based on these measures, a 
60 \eqn{\times} 60 correlation matrix was computed for each individual, subsequently converted into a dissimilarity matrix. In total, we have 20 such dissimilarity matrices. As additional external scales, NeuralScales gives 16 dimensions extracted from the psychological literature as a starting point for developing a theory of mental state representation: valence, arousal, warmth, competence, agency, experience, emotion, reason, mind, body, social, nonsocial, shared, and unique. 

}
\usage{
data("NeuralActivity")
data("NeuralScales")
data("NeuralScanner")
}
\format{
  A list of 20 dissimilarity matrices (NeuralActivity). 
  
  External scales (based on a questionnaire) containing proportions telling us to which degree people associate each of the 60 mental states to the 16 theoretical dimensions they extracted from the literature (NeuralScales).
  
  Scanner information on states, onset times and stimulus duration (NeuralScanner). 
  
  Head motion parameters (NeuralHM). 
}

\source{
Tamir D. I., Thornton M. A., Contreras J. M., & Mitchell J. P. (2015) Neural evidence that three dimensions organize mental state representation: rationality, social impact, and valence. Proceedings of the National Academy of Sciences of the United States of America, 113(1), 194-199.
}

\examples{
data(NeuralActivity)
str(NeuralActivity)

data(NeuralScales)
str(NeuralScales)

data(NeuralScanner)
str(NeuralScanner)
}
\keyword{datasets}
