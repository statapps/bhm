\name{print}
\alias{print.bhm}
\alias{print.brm}
\alias{print.mpl}
\alias{print.picreg}
\alias{print.pIndex}
\alias{print.resboot}
\alias{print.summary.bhm}
\alias{print.summary.brm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print a fitted object or a summary of fitted object
}
\description{
print and summary are used to provide a short summary of outputs from "bhm", "brm", "mpl", "pIndex", "resboot".
}
\usage{
\method{print}{bhm}(x, ...)
\method{print}{brm}(x, digits = 4, ...)
\method{print}{mpl}(x, digits = 3, ...)
\method{print}{pIndex}(x, ...)
\method{print}{picreg}(x, digits=3, ...)
\method{print}{resboot}(x, ...)
\method{print}{summary.bhm}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class returned from bhm, pIndex or resboot fit}
  \item{digits}{number if digits to be printed}
  \item{\dots}{other options used in print()}
}
\details{
  print.bhm is called to print object or summary of object from the biomarker threshold models \code{\link{bhm}}.  
  print.pIndex is called to print object or summary of object from the probability index model \code{\link{pIndex}}.
  print.resboot is called to print object or summary of object from the residuall bootstrap method for biomarker threshold models \code{\link{resboot}}.
  summary(fit) provides detail summary of `bhm' model fit, including parameter estimates, 
standard errors, and 95 percent CIs.

  The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}
\author{Bingshu E. Chen}
\seealso{
 The default method for print \code{\link{print.default}}. Other methods include
\code{\link{glm}}, 
\code{\link{bhm}},
\code{\link{brm}}, 
\code{\link{mpl}},
\code{\link{pIndex}},
\code{\link{resboot}}.
}
\examples{
#
#  print(fit)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ summary }
