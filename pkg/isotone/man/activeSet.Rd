\name{activeSet}
\alias{activeSet}
\alias{print.activeset}
\alias{summary.activeset}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Active Set Methods for Isotone Optimization}
\description{Isotone optimization can be formulated as a convex programming problem with simple linear constraints.
This functions offers active set strategies for a collection of isotone optimization problems pre-specified in the 
package. 
}
\usage{
activeSet(z, isomat, mySolver = lsSolver, ups = 1e-12, check = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing observed predictor values (see details for the response values)}
  \item{isomat}{Matrix with 2 columns that contains isotonicity conditions (see examples)}
  \item{mySolver}{Various functions are pre-defined (see details). For user-specified functions \code{fSolver} with additional 
  arguments can be used (see details as well).}
  \item{ups}{Upper boundary}
  \item{check}{If TRUE, KKT feasibility checks for isotonicity of the solution are performed}
  \item{...}{Additional arguments for \code{fSolver} (see details)}
}
\details{The observed response values \code{y} are passed to the solver. The following solvers are specified: \code{aSolver()} for asymmetric least squares, \code{dSolver()} for the least 
absolute value, \code{eSolver()} minimizes l1-approximation, \code{fSolver()} for arbitrary differentiable functions (additional
arguments \code{fobj} and \code{gobj} must be provided), \code{hSolver()} for Huber loss function, \code{iSolver()} for 
SILF loss (support vector regression), \code{lfSolver()} for least squares with non-diagonal weights, \code{lsSolver()} for least 
squares with diagonal weights (with response values \code{y} and weights \code{w} as additional arguments) , \code{mSolver()} for Chebyshev norm, \code{oSolver()} for power norm, \code{pSolver()} for quantile 
loss function, and finally \code{sSolver()} for Poisson likelihood. 

See corresponding help files and package vignette for additional description.
}

\value{
  Generates an object of class \code{activeset}.
  \item{z}{Vector containing the fitted values}
  \item{lambda}{Vector with Lagrange multipliers}
  \item{fval}{Value of the target function}
  \item{constr.vals}{Vector with the values of isotonicity constraints}
  \item{Alambda}{Constraint matrix multiplied by lambda (should be equal to gradient)}
  \item{gradient}{Gradient}
  \item{isocheck}{List containing the KKT checks for stationarity, primal feasibility, dual feasibility, and complementary slackness}
  \item{niter}{Number of iterations}
  \item{call}{Matched call}
}
\references{
de Leeuw, J., Hornik, K., Mair, P. (2008). Isotone optimization in R: Active Set methods and pool-adjacent-violators algorithm. 
Journal of Statistical Software, forthcoming.
}
\author{Jan de Leeuw, Kurt Hornik, Patrick Mair}

\seealso{\code{\link{gpava}}, \code{\link{lsSolver}}, \code{\link{dSolver}}, \code{\link{mSolver}}, \code{\link{fSolver}},
\code{\link{pSolver}}, \code{\link{lfSolver}}, \code{\link{oSolver}}, \code{\link{aSolver}}, \code{\link{eSolver}},
\code{\link{sSolver}}, \code{\link{hSolver}}, \code{\link{iSolver}}
}
\examples{

## Data specification
z <- 9:1                    ##predictor values (only order is relevant)
set.seed(12345)
y <- rnorm(9)               ##normal distributed response values
w1 <- rep(1,9)              ##unit weights
Atot <- cbind(1:8, 2:9)     ##Matrix defining isotonicity (total order)
Atot


## Least squares solver (pre-specified and user-specified)
fit.ls1 <- activeSet(z, Atot, lsSolver, y = y, weights = w1)
fit.ls1
summary(fit.ls1)
fit.ls2 <- activeSet(z, Atot, fSolver, fobj = function(x) sum(w1*(x-y)^2), gobj = function(x) 2*drop(w1*(x-y)))

## LS vs. GLS solver (needs weight matrix)
set.seed(12345)
wvec <- 1:9
wmat <- crossprod(matrix(rnorm(81),9,9))/9  
fit.wls <- activeSet(z, Atot, lsSolver, y = y, weights = wvec)
fit.gls <- activeSet(z, Atot, lfSolver, y = y, weights = wmat)


## Quantile regression
fit.qua <- activeSet(z, Atot, pSolver, y = y, weights = wvec, aw = 0.3, bw = 0.7)


## Mean absolute value norm
fit.abs <- activeSet(z, Atot, dSolver, y = y, weights = w1)

## Lp norm
fit.pow <- activeSet(z, Atot, oSolver, y = y, weights = w1, p = 1.2)

## Chebyshev norm
fit.che <- activeSet(z, Atot, mSolver, y = y, weights = w1)

## Efron's asymmetric LS
fit.asy <- activeSet(z, Atot, aSolver, y = y, weights = w1, aw = 2, bw = 1)

## Huber and SILF loss
fit.hub <- activeSet(z, Atot, hSolver, y = y, weights = w1, eps = 1)
fit.svm <- activeSet(z, Atot, iSolver, y = y, weights = w1, beta = 0.8, eps = 0.2)


## Negative Poisson log-likelihood
set.seed(12345)
yp <- rpois(9,5)
fit.poi <- activeSet(z, Atot, sSolver, y = yp)

## LS on tree ordering
Atree <- matrix(c(1,1,2,2,2,3,3,8,2,3,4,5,6,7,8,9),8,2)
Atree
fit.tree <- activeSet(z, Atree, lsSolver, y = y, weights = w1)


## LS on loop ordering
Aloop <- matrix(c(1,2,3,3,4,5,6,6,7,8,3,3,4,5,6,6,7,8,9,9),10,2)
Aloop
fit.loop <- activeSet(z, Aloop, lsSolver, y = y, weights = w1)


## LS on block ordering
Ablock <- cbind(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)),c(rep(c(4,5,6),3),rep(c(7,8,9),3)))
Ablock
fit.block <- activeSet(z, Ablock, lsSolver, y = y, weights = w1)

## Isotone LS regression using gpava and active set (same results)
pava.fitted <- gpava(z, y)$x
aset.fitted <- activeSet(z, Atot, lsSolver, weights = w1, y = y)$x
mse <- mean((pava.fitted - aset.fitted)^2)
mse
}
\keyword{models}
