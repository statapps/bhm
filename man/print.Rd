\name{print}
\alias{print.bhm}
\alias{summary.bhm}
\alias{print.summary.bhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print a bhm object or a summary of bhm object
}
\description{
print and summary are used to provide a short summary of bhm outputs.
}
\usage{
\method{print}{bhm}(x, ...)
\method{print}{summary.bhm}(x, ...)
\method{summary}{bhm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a bhm class returned from bhm fit}
  \item{\dots}{other options used in print()}
  \item{object}{object returned from bhm fit}
}
\details{
  print.bhm is called to print object or summary of object from the biomarker threshold model \code{\link{bhm}}.

  The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}
\author{Bingshu E. Chen and Fang Tian}
\seealso{
 The default method for print \code{\link{print.default}}.
\code{\link{glm}}
\code{\link{bhm}}
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
