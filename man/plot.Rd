\name{plot}
\alias{plot.bhm}
\alias{plot.brm}
\alias{plot.pIndex}
\alias{plot.resboot}
\alias{plot.residuals.brm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a fitted biomarker threhold model
}
\description{
Several different type of plots can be produced for biomarker threshold mdels. 
Plot method is used to provide a summary of outputs from "bhm", "pIndex", "resboot".

Use "methods(plot)" and the documentation for these for other plot methods.
}
\usage{
\method{plot}{bhm}(x, type = c("profile", "density"), ...)
\method{plot}{brm}(x, type = c("HR"), ...)
\method{plot}{pIndex}(x, ...)
\method{plot}{resboot}(x, ...)
\method{plot}{residuals.brm}(x, type="Martingale", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class returned from "bhm", "pIndex" or "resboot" fit.}
  \item{type}{type of plot in bhm object, "profile" to plot profile likelihood,   "density" to plot trace and density of the threshold distribution. "HR" to plot hazard ratio of the "brm" boject.}
  \item{\dots}{other options used in plot().}
}
\details{
  plot.bhm is called to plot either the profilelihood function or the threshold
  density function. 

  plot.pIndex is called to plot local probability index (\code{\link{pIndex}})
  of a continuous biomarker.

  plot.resboot is called to plot the bootstrap distribution of the likelihood 
  ratio test statistics for biomarker threshold models (\code{\link{resboot}}).

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
######## plot for bhm object
#
#  plot(fit, type = 'density')
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ summary }
