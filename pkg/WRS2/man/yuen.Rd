\name{yuen}
\alias{yuen}
\alias{yuenbt}
\alias{yuenbt}
\alias{pb2gen}
\alias{medpb2}
\alias{qcomhd}
\alias{print.yuen}
\alias{print.pb2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Independent samples t-tests on robust location measures.
}
\description{
The function \code{yuen} performs Yuen's test for trimmed means, \code{yuenbt} is a bootstrap version of it. The \code{pb2gen} function performs a t-test based on various robust estimators, \code{medpb2} compares two independent groups using medians, \code{qcomhd} compares arbitrary quantiles. 
}
\usage{
yuen(formula, data, tr = 0.2)
yuenbt(formula, data, tr = 0.2, nboot = 599, side = TRUE)
pb2gen(formula, data, est = "mom", nboot = 599)
medpb2(formula, data, nboot = 2000)
qcomhd(formula, data, q = c(0.1, 0.25, 0.5, 0.75, 0.9), nboot = 2000, alpha = 0.05, ADJ.CI = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula.
}
  \item{data}{
an optional data frame for the input data.
}
  \item{tr}{
trim level for the mean.
}
 \item{nboot}{
number of bootstrap samples.
}
 \item{side}{
 \code{side = TRUE} indicates two-sided method using absolute value of the  test statistics within the bootstrap; otherwise the equal-tailed method is used.
 }
\item{est}{
estimate to be used for the group comparisons: either \code{"onestep"} for one-step M-estimator of location using Huber's Psi, \code{"mom"} for the modified one-step (MOM) estimator of location based on Huber's Psi, or \code{"median"}, \code{"mean"}. 
}
\item{q}{
quantiles to be used for comparison.
}
\item{alpha}{
alpha level.
}
\item{ADJ.CI}{
whether CIs should be adjusted.
}
}

\details{
If \code{yuenbt} is used, p-value computed only when \code{side = TRUE}. \code{medpb2} is just a wrapper function for \code{pb2gen} with the median 
as M-estimator. It is the only known method to work well in simulations when tied values are likely to occur.\code{qcomhd} returns p-values and critical p-values based on Hochberg's method.
}

\value{
Returns objects of classes \code{"yuen"} or \code{"pb2"} containing:

\item{test}{value of the test statistic (t-statistic)}
\item{p.value}{p-value}
\item{conf.int}{confidence interval}
\item{df}{degress of freedom}
\item{diff}{trimmed mean difference}
\item{call}{function call}

\code{qcomhd} returns an object of class \code{"robtab"} containing:

\item{partable}{parameter table}
}
\references{
Wilcox, R. (2012). Introduction to Robust Estimation and Hypothesis Testing (3rd ed.). Elsevier.

Yuen, K. K. (1974). The two sample trimmed t for unequal population variances. Biometrika, 61, 165-170.
}
\seealso{
\code{\link{t1way}},\code{\link{t1waybt}}
}
\examples{
## Yuen's test
yuen(Anxiety ~ Group, data = spider)

## Bootstrap version of Yuen's test (symmetric CIs)
yuenbt(Anxiety ~ Group, data = spider)

## Using an M-estimator
pb2gen(Anxiety ~ Group, data = spider, est = "mom")
pb2gen(Anxiety ~ Group, data = spider, est = "mean")
pb2gen(Anxiety ~ Group, data = spider, est = "median")

## Using the median
medpb2(Anxiety ~ Group, data = spider)

## Quantiles
set.seed(123)
qcomhd(Anxiety ~ Group, data = spider, q = c(0.8, 0.85, 0.9), nboot = 500)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

