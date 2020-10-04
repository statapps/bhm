##########Note: git pull before working on the code.
\name{rpicexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{dpicexp}
\alias{ppicexp}
\alias{qpicexp}
\alias{rpicexp}

\alias{hpicexp}
\alias{Hpicexp}
\alias{picfit}
\alias{summary.picreg}

\title{
The Piecewise Exponential Distribution
}
\description{
Density, distribution function, quantile function, hazard function h(t), 
cumulative hazard function H(t), and random 
generation for the piecewise exponential distribution with rate equal to 'rate'
and cut points equal to 'cuts'.
}

\usage{
dpicexp(x, rate=1, cuts=c(0, 10), log = FALSE)
ppicexp(q, rate=1, cuts=c(0, 10), lower.tail = TRUE, index = NULL)
qpicexp(p, rate=1, cuts=c(0, 10), lower.tail = TRUE)
rpicexp(n, rate=1, cuts=c(0, 10))

hpicexp(x, rate, cuts, index=NULL)
Hpicexp(x, rate, cuts, index=NULL)
#
## to fit a piece exponential survival model use:
#
picfit = function(y, cuts=c(0, 10))
#
}

\arguments{
  \item{x, q}{vector of quantiles.}

  \item{p}{vector of probabilities.}

  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.}
  
  \item{rate}{vector rate parameter, defaulting to 1.}

  \item{cuts}{cut points, defaulting 0 to 10.}
  
  \item{log}{logical; if TRUE, probability p are given as log(p).}
  
  \item{lower.tail}{logical; if TRUE(default), probabilities are P[X <= x],
  otherwise, P[X>x].}

  \item{index}{index of x, q in the interval defined by cuts, it saves time if index is known.}
}

\details{
If the rate is not specified, it assumes the default value of 1.
}

\value{
%%  ~Describe the value returned
dpicexp gives the density, ppicexp gives the distribution function, 
qpicexp gives the quantile function, and 
rpicexp generates random deviates. 

The length of the result is determined by n for rpicexp.

Only the first elements of the logical arguments are used.
}

\references{
Chen, B. E., Cook, R. J., Lawless, J. F. and Zhan, M. (2005). 
Statistical methods for multivariate interval-censored recurrent events. 
Statistics in Medicine. Vol 24, 671-691.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}

\note{
The cumulative hazard H(t) = -log(1-F(t)) is log(1-ppicexp(t, rate, cuts)), 
or more efficiently call function Hpicexp(t, rate, cuts).
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{exp}} for the exponential function.

\code{\link{Distributions}} for other standard distributions, including 
\code{\link{dgamma}} for the gamma distribution and 
\code{\link{dweibull}} for the Weibull distribution.
}
\examples{
##
### No run
# n = 100
# rate = c(1, 1, 0.5, 0.25) 
# cuts = c(0, 1, 2.5, 5, 10)
# x = rpicexp(n, rate, cuts)
#
### compare rexp and rpicexp
#
#print(ppicexp(2.5, rate = .5))
#print(pexp(2.5, rate = 0.5))
#
#
}
\keyword{Piecewise exponential distribution}
\keyword{Piecewise exponential regression}

