\name{maxwell}
\alias{maxwell}
\docType{data}
\title{Maxwell's data}
\description{This data set is a hypothetical data set originally contrived by Maxwell (1961) 
for demonstrating his method of discriminant analysis. The data consist of three criterion groups, schizophrenic , manic-depressive
and anxiety state, and four binary predictor variables each indicating either
presence (1) or absence (0) of a certain symptom. The four symptoms are anxiety
suspicion, schizophrenic type of thought disorders, and delusions of guilt.
These four binary variables were factorially combined to form 16 distinct patterns of symptoms
(predictor patterns), and each of these patterns is identified with a row of the table, which
contains the cross-classification of 620 patients according to the 16 patterns of symptoms
and the three criterion groups.
}
\usage{data(maxwell)}
\format{A list with the frequency table as the first element and the row covariates as the second.
}
\details{This dataset can be used for canonical CA. The binary predictor variables can
be considered as row covariates.
}

\references{Maxwell, A.E. (1961). Canonical variate analysis when the variables are dichotomous. 
Educational and Psychological Measurement, 21,259-271.
}
\examples{
data(maxwell)
## maybe str(maxwell) ; plot(maxwell) ...
}
\keyword{datasets}
