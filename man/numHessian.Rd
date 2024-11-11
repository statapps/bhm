\name{numHessian}
\alias{numHessian}

\title{Calculate Hessian or Information Matrix}
\description{Calculate a numerical approximation to the Hessian matrix of a 
function at a parameter value.}
\usage{
    numHessian(func, theta, h = 0.0001, method=c("fast", "easy"), ...)
}
\arguments{
    \item{func}{a function for which the first (vector) argument 
        is used as a parameter vector.}
    \item{theta}{the parameter vector first argument to func.}
    \item{h}{the step used in the numerical calculation.}
    \item{method}{one of "fast" or "easy" indicating 
    the method to use for the approximation.}
    \item{...}{additional named or unmaned arguments to be passed to \code{func}.}
}
\value{An p by p matrix of the Hessian of the function calculated at the 
       point \code{theta}. If the \code{func} is a log likelihood function, 
       then the negative of the p by p matrix is the information matrix.}

\details{
  The function \code{numHessian} calculates an numerical approximation to
  the p by p second order derivative of a scalar real valued function with p-vector
  argument theta.
  This function can be used to check if the information matrix of a log likelihood is correct or not.
}

\seealso{
    \code{\link{numScore}}
}

\examples{
  g = function(x, a) (x[1]+2*x[2]^3 - x[3]^3 + a*sin(x[1]*x[2]))
  x0= c(1, 2, 3)
  numHessian(g, theta = x0, a = 9)
  numHessian(g, theta = x0, method = 'easy', a = 9)
}
