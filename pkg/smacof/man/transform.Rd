\name{transform}
\alias{transform}

\title{Internal Dissimilarity Transformation}
\description{Function that performs the optimal scaling calls (used internally)
}
\usage{
transform(Target, x, w = rep(1,length(x$x)), normq = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{Target}{unconstrained vector of target values}
   \item{x}{object of type optScal}
   \item{w}{vector non-negative weights}
   \item{normq}{sum of squares normalization}
}  
