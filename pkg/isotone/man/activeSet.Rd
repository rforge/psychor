\name{activeSet}
\alias{activeSet}
\alias{print.activeset}
\alias{aSolver}
\alias{dSolver}
\alias{eSolver}
\alias{fSolver}
\alias{hSolver}
\alias{iSolver}
\alias{lfSolver}
\alias{lsSolver}
\alias{mSolver}
\alias{oSolver}
\alias{pSolver}
\alias{sSolver}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Active Set Methods for Isotone Optimization}
\description{Isotone optimization can be formulated as a convex programming problem with simple linear constraints.
This functions offers active set strategies for a collection of isotone optimization problems pre-specified in the 
package. 
}
\usage{
activeSet(x, isomat, mySolver = lsSolver, ups = 1e-12, check = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector containing predictors}
  \item{isomat}{Matrix with 2 columns that contains isotonicity conditions (see examples)}
  \item{mySolver}{Various functions are pre-defined (see details). For user-specified functions \code{fSolver} with additional 
  arguments can be used (see details as well).}
  \item{ups}{Upper boundary}
  \item{check}{If TRUE, KKT feasibility checks for isotonicity of the solution are performed}
  \item{...}{Additional arguments for \code{fSolver} (see details)}
}
\details{The following solvers are specified: \code{aSolver()} for asymmetric least squares, \code{dSolver()} for the least 
absolute value, \code{eSolver()} minimizes l1-approximation, \code{fSolver()} for arbitrary differentiable functions (additional
arguments \code{fobj} and \code{gobj} must be provided; see examples), \code{hSolver()} for Huber loss function, \code{iSolver()} for 
SILF loss (support vector regression), \code{lfSolver()} for least squares with non-diagonal weights, \code{lsSolver()} for least 
squares with diagonal weights (with response values \code{z} and weights \code{w} as additional arguments) , \code{mSolver()} for Chebyshev norm, \code{oSolver()} for power norm, \code{pSolver()} for quantile 
loss function, and finally \code{sSolver()} for Poisson likelihood. 

More descriptions can be found in the package vignette.  
}

\value{
  Generates an object of class \code{activeset}.
  \item{x}{Vector containing the fitted values}
  \item{lambda}{Vector with Lagrange multipliers}
  \item{func.vals}{Value of the target function}
  \item{constr.vals}{Vector with the values of isotonicity constraints}
  \item{Alambda}{Constraint matrix multiplied by lambda (should be equal to gradient)}
  \item{gradient}{Gradient}
  \item{isocheck}{List containing the KKT checks for stationarity, primal feasibility, dual feasibility, and complementary slackness}
  \item{call}{Matched call}
}
\references{
de Leeuw, J., Hornik, K., Mair, P. (2008). Isotone optimization in R: Active Set methods and pool-adjacent-violators algorithm. 
Journal of Statistical Software, forthcoming.
}
\author{Jan de Leeuw, Kurt Hornik, Patrick Mair}


\examples{

##Example 1: Comparison with PAVA
set.seed(12345)
x0 <- 9:1                   ##predictor values
z <- rnorm(9)               ##response values
wu <- rep(1,9)              ##unit weights
btota <- cbind(1:8, 2:9)    ##Matrix defining isotonicity (total order)

##Now we apply PAVA, LS active set, and, again LS active set by specifying function (fobj) and derivative (gobj)          
fit.pava <- gpava(x0, z)
fit.ls <- activeSet(x0, btota, lsSolver, w = wu, z = z)
fit.f <- activeSet(x0, btota, fSolver, fobj = function(x) sum(wu*(x-z)^2), gobj = function(x) 2*drop(wu*(x-z)))

##Matrix with fitted values (basically same results for each method)
val.mat <- cbind(fit.pava$yfit, fit.ls$x, fit.f$x)
val.mat


##Example 2: Diagonal vs. non-diagonal weights
set.seed(12345)
ww <- 1:9                                  ##weight vector
wf <- crossprod(matrix(rnorm(81),9,9))/9   ##symmetric weight matrix 

activeSet(x0, btota, lsSolver, check = TRUE, w = ww, z = z)

#FIXME!!! activeSet(x0, btota, lfSolver, check = TRUE, x0, w = wf, z = z) 

##Example 3: Other norms
##Approximate l1 with additional argument eps = 1e-4
activeSet(x0, btota, eSolver, check = TRUE, z = z, w = wu, eps = 1e-4)

##Approximate l1 with additional power argument p = 1.2
activeSet(x0, btota, oSolver, check = TRUE, z = z, w = wu, p = 1.2)

##Least absolute values (exact l1)
activeSet(x0, btota, dSolver, check = TRUE, w = wu, z = z)

##Approximate l(infty) with additional power argument p = 7
#FIXME!!! activeSet(x0, btota, oSolver, check = TRUE, z = z, w = wu, p = 7)

##Exact l(infty) (Chebyshew norm)
activeSet(x0, btota, mSolver, check = TRUE, w = wu, z = z)

##Poisson likelihood
set.seed(12345)    
z <- rpois(9, 5)
activeSet(x0, btota, sSolver, check = TRUE, w = wu, z = z)

##Asymmetric least squares with additional asymmetric weights aw and bw
activeSet(x0, btota, aSolver, check = TRUE, z=z, w = wu, aw = 2, bw = 1)

##Huber norm with additional argument eps = 1
activeSet(x0, btota , hSolver, check = TRUE, z = z, w = wu, eps = 1)

##SILF norm with additional arguments beta = 0.8 and eps = 0.2
activeSet(x0, btota, iSolver, check = TRUE, z = z, w = wu, beta = 0.8, eps = 0.2)


##Example 4: Other orders
## Tree order
btree <- matrix(c(1,1,2,2,2,3,3,8,2,3,4,5,6,7,8,9),8,2)
activeSet(x0, btree, lsSolver, check = TRUE, w = wu, z = z)

## Block order 
#bprim <- cbind(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)),
#             c(rep(c(4,5,6),3),rep(c(7,8,9),3)))
#activeSet(x0, bprim, lsSolver, check = TRUE, w = wu, z = z)

## Loop order
#bloop <- matrix(c(1,2,3,3,4,5,6,6,7,8,3,3,4,5,6,6,7,8,9,9),10,2)
#activeSet(x0, bloop, lsSolver, check = TRUE, w = wu, z = z)

}
\keyword{models}
