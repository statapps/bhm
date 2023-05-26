\name{multiRoot}
\alias{multiRoot}

\title{m-Dimensional Root (Zero) Finding}
\description{The function \code{multiRoot} searches for roots (i.e, zeros) of 
the vector-valued function \code{func} with repsect to its first argument.
} 
\usage{
    multiRoot(func, theta, ..., verbose = FALSE, maxIter = 20, tol = .Machine$double.eps^0.25)
}
\arguments{
    \item{func}{a m-vector function for which the roots is sought.}
    \item{theta}{the parameter vector first argument to func.}
    \item{...}{an additional named or unmaned arguments to be passed to \code{func}.}
    \item{verbose}{print out the verbose, default is FALSE. }
    \item{maxIter}{the maximum number of iterations, default is 20.}
    \item{tol}{the desired accuracy (convergence tolerance), default is .Machine$double.eps^0.25.}
}
\value{A list with at least four components:
  
  \item{root}{a vector of theta that solves func(theta) = 0.}
  \item{r.root}{a vector of f(roots) that evaluates at theta = roots.}
  \item{iter}{number of iteratins used in the algorithm.}
  \item{convergence}{1 if the algorithm converges, 0 otherwise.}
}

\details{
  The function \code{multiRoot} finds an numerical approximation to
  func(theta) = 0. 
  This function can be used to solve the score function for a maximum 
  log likelihood estimate.

  This function make use of \code{numJacobian} calculates an numerical approximation to
  the m by p first order derivative of a m-vector real valued function.
  The parameter theta is updated by the Newton-Ralphson method

  theta = theta - solve((t(J) x J), J x func(theta))

  When m > p,  if the nonlinear system has not solution, 
  the method attempts to find a solution in the non-linear least squares sense.
  The sum of square sum(t(U)xU), where U = func(theta), will be minized. 
}

\seealso{
    \code{\link{numJacobian}}
    \code{\link{numScore}}
}

\examples{
  g = function(x, a) (c(x[1]+2*x[2]^3, x[2] - x[3]^3, a*sin(x[1]*x[2])))
  theta = c(1, 2, 3)
  multiRoot(g, theta, a = -3)
}
