\name{plot}
\alias{plot.pIndex}
\alias{plot.resboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot a fitted object or a summary of fitted object
}
\description{
plot and summary are used to provide a short summary of outputs from "bhm", "pIndex", "resboot".
}
\usage{
\method{plot}{pIndex}(x, ...)
\method{plot}{resboot}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class returned from "pIndex" or "resboot" fit}
  \item{\dots}{other options used in plot()}
}
\details{
  plot.pIndex is called to plot object or summary of object from the probability index model \code{\link{pIndex}}.
  plot.resboot is called to plot object or summary of object from the residuall bootstrap method for biomarker threshold models \code{\link{resboot}}.
  summary(fit) provides detail summary of `bhm' model fit, including parameter estimates, 
standard errors, and 95 percent CIs.

  The default method, plot.default has its own help page. Use methods("plot") to get all the methods for the plot generic.
}
\author{Bingshu E. Chen}
\seealso{
 The default method for plot \code{\link{plot.default}}.
\code{\link{glm}}
\code{\link{bhm}}
\code{\link{pIndex}}
\code{\link{resboot}}
}
\examples{
#
#  plot(fit)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ summary }
