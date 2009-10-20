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
activeSet(isomat, mySolver = "LS", x0 = NULL, ups = 1e-12, check = TRUE, maxiter = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isomat}{Matrix with 2 columns that contains isotonicity conditions, i.e. for row i it holds that fitted value i column 1 <= fitted value i column 2 (see examples)}
  \item{mySolver}{Various functions are pre-defined (see details). Either to funtction name or the corresponding string equivalent can be used. For user-specified functions \code{fSolver} with additional 
  arguments can be used (see details as well).}
  \item{x0}{Feasible starting solution. If \code{NULL} the null-vector is used internally.}
  \item{ups}{Upper boundary}
  \item{check}{If TRUE, KKT feasibility checks for isotonicity of the solution are performed}
  \item{maxiter}{Iteration limit}
  \item{...}{Additional arguments for the various solvers (see details)}
}
\details{The following solvers are specified. Note that \code{y} as the vector of observed values and \code{weights} as the vector of weights need to provided through \code{...} for each solver (except for \code{fSolver()} and \code{sSolver()}). Some solvers need additional arguments as described in the corresponding solver help files. More technical details can be found in the package vignette.    

The pre-specified solvers are the following (we always give the corresponding string equivalent in brackets): 
\code{lsSolver()} (\code{"LS"}) for least squares with diagonal weights, \code{aSolver()} (\code{"asyLS"}) for asymmetric least squares, \code{dSolver()} (\code{"L1"}) for the least absolute value, \code{eSolver()} (\code{"L1eps"}) minimizes l1-approximation. \code{hSolver()} (\code{"huber"}) for Huber loss function, \code{iSolver()} (\code{"SILF"}) for SILF loss (support vector regression), \code{lfSolver()} (\code{"GLS"}) for general least squares with non-diagonal weights, \code{mSolver()} (\code{"chebyshev"}) for Chebyshev L-inf norm, \code{oSolver()} (\code{"Lp"}) for L-p power norm, \code{pSolver()} (\code{"quantile"}) for quantile loss function, and finally \code{sSolver()} (\code{"poisson"}) for Poisson likelihood. 

\code{fSolver()} for user-specified arbitrary differentiable functions. The arguments \code{fobj} (target function ) and \code{gobj} (first derivative) must be provided plus any additional arguments used in the definition of \code{fobj}. 

}

\value{
  Generates an object of class \code{activeset}.
  \item{x}{Vector containing the fitted values}
  \item{y}{Vector containing the observed values}
  \item{lambda}{Vector with Lagrange multipliers}
  \item{fval}{Value of the target function}
  \item{constr.vals}{Vector with the values of isotonicity constraints}
  \item{Alambda}{Constraint matrix multiplied by lambda (should be equal to gradient)}
  \item{gradient}{Gradient}
  \item{isocheck}{List containing the KKT checks for stationarity, primal feasibility, dual feasibility, and complementary slackness (>= 0 means feasible)}
  \item{niter}{Number of iterations}
  \item{call}{Matched call}
}
\references{
de Leeuw, J., Hornik, K., Mair, P. (2009). Isotone optimization in R: Active Set methods and pool-adjacent-violators algorithm. 
Journal of Statistical Software, 32(5), 1-24.
}
\author{Jan de Leeuw, Kurt Hornik, Patrick Mair}

\seealso{\code{\link{gpava}}, \code{\link{lsSolver}}, \code{\link{dSolver}}, \code{\link{mSolver}}, \code{\link{fSolver}},
\code{\link{pSolver}}, \code{\link{lfSolver}}, \code{\link{oSolver}}, \code{\link{aSolver}}, \code{\link{eSolver}},
\code{\link{sSolver}}, \code{\link{hSolver}}, \code{\link{iSolver}}
}
\examples{

## Data specification
set.seed(12345)
y <- rnorm(9)               ##normal distributed response values
w1 <- rep(1,9)              ##unit weights
Atot <- cbind(1:8, 2:9)     ##Matrix defining isotonicity (total order)
Atot


## Least squares solver (pre-specified and user-specified)
fit.ls1 <- activeSet(Atot, "LS", y = y, weights = w1)
fit.ls1
summary(fit.ls1)
fit.ls2 <- activeSet(Atot, fSolver, fobj = function(x) sum(w1*(x-y)^2), gobj = function(x) 2*drop(w1*(x-y)), y = y, weights = w1)

## LS vs. GLS solver (needs weight matrix)
set.seed(12345)
wvec <- 1:9
wmat <- crossprod(matrix(rnorm(81),9,9))/9  
fit.wls <- activeSet(Atot, "LS", y = y, weights = wvec)
fit.gls <- activeSet(Atot, "GLS", y = y, weights = wmat)


## Quantile regression
fit.qua <- activeSet(Atot, "quantile", y = y, weights = wvec, aw = 0.3, bw = 0.7)


## Mean absolute value norm
fit.abs <- activeSet(Atot, "L1", y = y, weights = w1)

## Lp norm
fit.pow <- activeSet(Atot, "Lp", y = y, weights = w1, p = 1.2)

## Chebyshev norm
fit.che <- activeSet(Atot, "chebyshev", y = y, weights = w1)

## Efron's asymmetric LS
fit.asy <- activeSet(Atot, "asyLS", y = y, weights = w1, aw = 2, bw = 1)

## Huber and SILF loss
fit.hub <- activeSet(Atot, "huber", y = y, weights = w1, eps = 1)
fit.svm <- activeSet(Atot, "SILF", y = y, weights = w1, beta = 0.8, eps = 0.2)


## Negative Poisson log-likelihood
set.seed(12345)
yp <- rpois(9,5)
x0 <- 1:9
fit.poi <- activeSet(Atot, "poisson", x0 = x0, y = yp)

## LS on tree ordering
Atree <- matrix(c(1,1,2,2,2,3,3,8,2,3,4,5,6,7,8,9),8,2)
Atree
fit.tree <- activeSet(Atree, "LS", y = y, weights = w1)


## LS on loop ordering
Aloop <- matrix(c(1,2,3,3,4,5,6,6,7,8,3,3,4,5,6,6,7,8,9,9),10,2)
Aloop
fit.loop <- activeSet(Aloop, "LS", y = y, weights = w1)


## LS on block ordering
Ablock <- cbind(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)),c(rep(c(4,5,6),3),rep(c(7,8,9),3)))
Ablock
fit.block <- activeSet(Ablock, "LS", y = y, weights = w1)

## Isotone LS regression using gpava and active set (same results)
pava.fitted <- gpava(y = y)$x
aset.fitted <- activeSet(Atot, "LS", weights = w1, y = y)$x
mse <- mean((pava.fitted - aset.fitted)^2)
mse
}
\keyword{models}
