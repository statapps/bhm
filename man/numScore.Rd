\name{numScore}
\alias{numScore}

\title{Calculate the Score Function}
\description{Calculate a numerical approximation to the Score function of a 
function at a parameter value.}
\usage{
    numScore(func, theta, h = 0.0001, ...)
}
\arguments{
    \item{func}{a function for which the first (vector) argument 
        is used as a parameter vector.}
    \item{theta}{the parameter vector first argument to func.}
    \item{h}{the step used in the numerical calculation.}
    \item{...}{an additional arguments passed to \code{func}.}
}
\value{An p by 1 vector of the score of the function calculated at the 
       point \code{theta}. If the \code{func} is a log likelihood function, 
       then the p by 1 vector is the score function.}

\details{
  The function \code{numScore} calculates an numerical approximation to
  the p by 1 first order derivative of a scalar real valued function with p-vector
  argument theta. 
  This function can be used to check if the score function of a log likelihood is correct or not.
}

\seealso{
    \code{\link{numHessian}}
}

\examples{
  g = function(x, a) (x[1]+2*x[2]^3 - x[3]^3 + a*sin(x[1]*x[2]))
  x0 = c(1, 2, 3)
  numScore(g, x0, a = -3)
}
