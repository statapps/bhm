\name{print}
\alias{print.bhm}
\alias{summary.bhm}
\alias{print.summary.bhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print an bhm object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Bingshu E. Chen and Fang Tian}
\note{
%%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{coxph}}
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
