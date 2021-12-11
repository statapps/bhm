\name{llm}
\alias{llm}
\alias{llm.default}
\alias{llm.formula}
\alias{plot.llm}
\alias{print.llm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a L-shape linear model}
\description{
Fit a L-shape linear model with the cut point estimated by the 
profie likelihood method.
}
\usage{
\method{llm}{formula}(formula, data=list(...), epsilon = 0.025, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "\code{formula}": a sympolic description of the model to be fitted. The details of model specificatoin are given under 'Details'.}
  \item{data}{an optional data frame, list or enviroment containing the variables in the model.}
  \item{epsilon}{Step width for the profile likelihood method, default is 0.025}
  \item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below).}
}
\details{
Define a L shape linear funcation that change slope at c0:

when x <c0, y = b1 + b2*x

when x>=c0, y = b1 + b2*x  + b3*(x-x0)
             = (b1-b3*c0) + (b2+b3)*x

}
\value{
llm returns an object of class inheriting from "llm" which inherits from the class glm. See later in this section.

An object of class "llm" is a list containing at least the following components:

  \item{coefficients}{a named vector of coefficients from 'llm'}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted.values}{the fitted mean values.}
  \item{rank}{the numeric rank of the fitted linear model.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{call}{the matched call.}
  \item{terms}{the 'terms' object used.}
  \item{c.max}{the maximum likelihood estimate for the threshold parameter(s).}
  \item{loglik}{the log-likelihood with the final values of the coefficients.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}
%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{brm}}, 
\code{\link{lm}},
\code{\link{glm}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

n = 20 ; p <- 2
X = matrix(rnorm(n * p), n, p) # no intercept!
w = X[, 1]; age = X[, 2]

wc = w - 0.2; sigma = 0.5
y = rnorm(n, -0.1+0.7*w-1.2*ifelse(wc>0, wc, 0), sigma)

fit=llm(y~w+age, epsilon = 0.005)
print(fit)
print(summary(fit))
#### No run
## plot(fit)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{linear model}
\keyword{threshold models}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
