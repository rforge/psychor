\name{activeSet}
\alias{activeSet}
\alias{print.activeset}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Active Set Methods for Isotone Optimization}
\description{Isotone optimization can be formulated as a convex programming problem with simple linear constraints.
This functions offers active set strategies for a collection of isotone optimization problems pre-specified in the 
package. 
}
\usage{
activeSet(z, isomat, mySolver = lsSolver, ups = 1e-12, check = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Vector containing predictors}
  \item{isomat}{Matrix with 2 columns that contains isotonicity conditions (see examples)}
  \item{mySolver}{Various functions are pre-defined (see details). For user-specified functions \code{fSolver} with additional 
  arguments can be used (see details as well).}
  \item{ups}{Upper boundary}
  \item{check}{If TRUE, KKT feasibility checks for isotonicity of the solution are performed}
  \item{...}{Additional arguments for \code{fSolver} (see details)}
}
\details{The following solvers are specified: \code{aSolver()} for asymmetric least squares, \code{dSolver()} for the least 
absolute value, \code{eSolver()} minimizes l1-approximation, \code{fSolver()} for arbitrary differentiable functions (additional
arguments \code{fobj} and \code{gobj} must be provided), \code{hSolver()} for Huber loss function, \code{iSolver()} for 
SILF loss (support vector regression), \code{lfSolver()} for least squares with non-diagonal weights, \code{lsSolver()} for least 
squares with diagonal weights (with response values \code{z} and weights \code{w} as additional arguments) , \code{mSolver()} for Chebyshev norm, \code{oSolver()} for power norm, \code{pSolver()} for quantile 
loss function, and finally \code{sSolver()} for Poisson likelihood. 

See corresponding help files and package vignette for additional description.
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

\seealso{\code{\link{gpava}}, \code{\link{lsSolver}}, \code{\link{dSolver}}, \code{\link{mSolver}}, \code{\link{fSolver}},
\code{\link{pSolver}}, \code{\link{lfSolver}}, \code{\link{oSolver}}, \code{\link{aSolver}}, \code{\link{eSolver}},
\code{\link{sSolver}}, \code{\link{hSolver}}, \code{\link{iSolver}}
}
\examples{

##Example 1: Comparison with PAVA
set.seed(12345)
z <- 9:1                   ##predictor values
y <- rnorm(9)               ##response values
w <- rep(1,9)              ##unit weights
btota <- cbind(1:8, 2:9)    ##Matrix defining isotonicity (total order)

##Now we apply PAVA, LS active set, and, again LS active set by specifying function (fobj) and derivative (gobj)          
fit.pava <- gpava(z, y)
#fit.ls <- activeSet(z, btota, lsSolver, weights = w, y = y)
fit.f <- activeSet(z, btota, fSolver, fobj = function(x) sum(w*(x-y)^2), gobj = function(x) 2*drop(w*(x-y)))

##Matrix with fitted values (basically same results for each method)
#val.mat <- cbind(fit.pava$yfit, fit.ls$x, fit.f$x)
#val.mat


##Example 2: Diagonal vs. non-diagonal weights
set.seed(12345)
ww <- 1:9                                  ##weight vector
wf <- crossprod(matrix(rnorm(81),9,9))/9   ##symmetric weight matrix 

#activeSet(z, btota, lsSolver, check = TRUE, weights = ww, y = y)

#FIXME!!! activeSet(z, btota, lfSolver, check = TRUE, z, weights = wf, y = y) 

##Example 3: Other norms
##Approximate l1 with additional argument eps = 1e-4
activeSet(z, btota, eSolver, check = TRUE, y = y, weights = w, eps = 1e-4)

##Approximate l1 with additional power argument p = 1.2
activeSet(z, btota, oSolver, check = TRUE, y = y, weights = w, p = 1.2)

##Least absolute values (exact l1)
#activeSet(z, btota, dSolver, check = TRUE, weights = w, y = y)

##Approximate l(infty) with additional power argument p = 7
#FIXME!!! activeSet(z, btota, oSolver, check = TRUE, y = y, weights = wu, p = 7)

##Exact l(infty) (Chebyshew norm)
#activeSet(z, btota, mSolver, check = TRUE, weights = w, y = y)

##Poisson likelihood
set.seed(12345)    
y <- rpois(9, 5)
#activeSet(z, btota, sSolver, check = TRUE, weights = w, y = y)

##Asymmetric least squares with additional asymmetric weights aw and bw
activeSet(z, btota, aSolver, check = TRUE, y=y, weights = w, aw = 2, bw = 1)

##Huber norm with additional argument eps = 1
activeSet(z, btota , hSolver, check = TRUE, y = y, weights = w, eps = 1)

##SILF norm with additional arguments beta = 0.8 and eps = 0.2
activeSet(z, btota, iSolver, check = TRUE, y = y, weights = w, beta = 0.8, eps = 0.2)


##Example 4: Other orders
## Tree order
btree <- matrix(c(1,1,2,2,2,3,3,8,2,3,4,5,6,7,8,9),8,2)
#activeSet(z, btree, lsSolver, check = TRUE, weights = w, y = y)

## Block order 
#bprim <- cbind(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)),
#             c(rep(c(4,5,6),3),rep(c(7,8,9),3)))
#activeSet(z, bprim, lsSolver, check = TRUE, weights = wu, y = y)

## Loop order
#bloop <- matrix(c(1,2,3,3,4,5,6,6,7,8,3,3,4,5,6,6,7,8,9,9),10,2)
#activeSet(z, bloop, lsSolver, check = TRUE, weights = wu, y = y)

}
\keyword{models}
