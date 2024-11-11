\name{multiRoot}
\alias{multiRoot}

\title{m-Dimensional Root (Zero) Finding}
\description{The function \code{multiRoot} searches for root (i.e, zero) of 
the vector-valued function \code{func} with respect to its first argument using 
the Gauss-Newton algorithm.
} 
\usage{
    multiRoot(func, theta, ..., verbose = FALSE, maxIter = 50, 
           thetaUp = NULL, thetaLow = NULL, tol = .Machine$double.eps^0.25)
}
\arguments{
    \item{func}{a m-vector function for which the root is sought.}
    \item{theta}{the parameter vector first argument to func.}
    \item{thetaLow}{the lower bound of theta.}
    \item{thetaUp}{the upper bound of theta.}
    \item{verbose}{print out the verbose, default is FALSE. }
    \item{maxIter}{the maximum number of iterations, default is 20.}
    \item{tol}{the desired accuracy (convergence tolerance), default is .Machine$double.eps^0.25.}
    \item{...}{an additional named or unmaned arguments to be passed to \code{func}.}
}

\value{A list with at least four components:
  
  \item{root}{a vector of theta that solves func(theta) = 0.}
  \item{f.root}{a vector of f(root) that evaluates at theta = root.}
  \item{iter}{number of iteratins used in the algorithm.}
  \item{convergence}{1 if the algorithm converges, 0 otherwise.}
}

\details{
  The function \code{multiRoot} finds an numerical approximation to
  func(theta) = 0 using Newton method: theta = theta - solve(J, func(theta)) when m = p. 
  This function can be used to solve the score function euqations for a maximum 
  log likelihood estimate.

  This function make use of \code{numJacobian} calculates an numerical approximation to
  the m by p first order derivative of a m-vector valued function.
  The parameter theta is updated by the Gauss-Newton method:

  theta = theta - solve((t(J) x J), J x func(theta))

  When m > p,  if the nonlinear system has not solution, 
  the method attempts to find a solution in the non-linear least squares sense 
  (Gauss-Newton algorithm).
  The sum of square sum(t(U)xU), where U = func(theta), will be minimized. 
}

\seealso{
    \code{\link{optim}} (which is preferred) and \code{\link{nlm}}, 
    \code{\link{nlminb}},
    \code{\link{numJacobian}}, 
    \code{\link{numScore}}, 
    \code{\link{optimize}} and \code{\link{uniroot}} for one-dimension optimization.
}

\references{
Gauss, Carl Friedrich(1809). Theoria motus corporum coelestium in sectionibus conicis solem ambientum.}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}

\examples{
  g = function(x, a) (c(x[1]+2*x[2]^3, x[2] - x[3]^3, a*sin(x[1]*x[2])))
  theta = c(1, 2, 3)
  multiRoot(g, theta, a = -3)
}
