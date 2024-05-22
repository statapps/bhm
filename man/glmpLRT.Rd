\name{glmpLRT}
\alias{glmpLRT}
\alias{glmpLRT.default}
\alias{glmpLRT.formula}
\alias{plot.glmpLRT}
\alias{print.glmpLRT}
\alias{summary.glmpLRT}
\title{Penalized likelihood ratio test for the generalized linear models.
}
\description{
\{glmpLRT\} is an R function for the penalized likelihood ratio test in generalized lienar models.
It uses the penalized likehood method (Gavanji, Jiang and Chen, 2021)
to identify a cut-point (thershold parameter) for the biomarker in the generalized linear models.
The model is specified by giving a symbolic description of the linear predictor and 
a description of the distribution family.
}
\usage{
%glmpLRT(x, ...)
\method{glmpLRT}{formula}(formula, family = binomial, data=list(...), lambda = 15, 
  c0 = 0.5, p1 = 1, method = c("pLRT", "Davies", "Bootstrap"), B=10, K = 50,
  epsilon = 0.025,...)
# use 
#          glmpLRT(y ~ biomarker)             
# or
#          glmpLRT(y ~ biomarker + x1 + x2, p1=0)
#
# to fit a prognostic model with biomarker term adjust for 
# covariates x1, x2, etc.
#
# use 
#
#          glmpLRT(y ~ biomarker+x1+x2+x3+x4+x5, p1=2...)
#
# to fit a predictive model with interaciton between biomarker 
# and x1, x2, adjusted for x3, x4, x5, etc.
# 
# use 
#          glmpLRT(x, y, control, ...)
#
# to fit a model without formula
#
# Biomarker shall be in the first dependent variable
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "\code{formula}": a sympolic description of the model to be fitted. The details of model specificatoin are given under 'Details'.}
  \item{family}{a description of the error distribution and the link function to be used in the glm model. (See \code{family} for details of family functions.)}
  \item{data}{an optional data frame, list or enviroment containing the variables in the model.}
  \item{method}{the method to be used to fit a glmpLRT model. The default method is the "pLRT" method for the penalized likelihood ratio test. Other options include "Davis" for the Davis test and "Bootstrap" for the bootstrap test.}
  \item{lambda}{the tuning parameter for the penalized likelihood function, default is lambda = 15.}
  \item{c0}{the cut point for the threshold function, default is c0 = 0.5. Other options include c0 = 0.25 and c0 = 0.75.}
  \item{p1}{the number of covariates interact with the biomarker variable, the default is p1 = 1. If p1 = 0, then fit a model with biomarker main effect, adjusted for other potential covariates.}
  \item{B}{the number of bootstrap sample if the Bootstrap method is used.}
  \item{K}{smoothing parameter for the indicator function, the large value of K the better approximation of the indicator function, default is K = 50.}
  \item{epsilon}{Step width for the profile likelihood method, default is epsilon = 0.025.}
  \item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below).}  

%  \item{x}{For "glmpLRT.default", x is a design matrix of dimension n * p.}
%  \item{y}{y is a vector of observations of length n for the glm models}
}
\details{
'biomarker' is a Biomarker variable. This variable is required and shall be the first dependent variable in the formula. It is not necessary to use the variable name 'biomarker', other variable name is allowed too.

'p1' controls the number of other covariates that interact with the biomarker variable. 

"glmpLRT.default" is the workhorse functions: they are not normally
called directly but can be more efficient where the response 
vector, design matrix and family have already been calculated.

"print", "plot" and "summary" methods can be applied to a fitted "glmpLRT" class to display the results.
}
\value{
glmpLRT returns an object of class inheriting from "glmpLRT" which inherits from the class glm. See later in this section.  

An object of class "glmpLRT" is a list containing at least the following components:
  
  \item{coefficients}{a named vector of coefficients from 'glmpLRT'}
  \item{c.max}{the maximum likelihood estimate for the threshold parameter(s).}
  \item{loglik}{the log-likelihood with the final values of the coefficients.}
  \item{linear.predictors}{the vector of linear predictors, one per subject.}
  \item{mpv}{p-value for the penalized likelihood ratio test.}
  \item{rpv}{p-value for the Davis test.}
  \item{bpv}{p-value for the Bootstrap test.}
}
\references{
Gavanji, P., Jiang, W. and Chen, B. E. (2021).
Penalized likelihood ratio test for a biomarker threshold effect in clinical trials based on generalized linear models.
Canadian Journal of Statistics. 
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}
%\note{
%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{glmpLRT}}, 
\code{\link{glm}},
\code{\link{plot.glmpLRT}},
\code{\link{print.glmpLRT}}
}
\examples{
## A simulated example
bmk = rnorm(100, 3, 0.25)
age = rnorm(100, 50, 20)
trt = rbinom(100, 1, 0.5)
lp = exp(log(0.25) + 0.1*ifelse(bmk>2.5, 1, 0) + 0.69*trt)
p = lp/(1+lp)
y = rbinom(100, 1, p)
fit = glmpLRT(y~bmk+trt+age, p1 = 0)
print(fit)
}
\keyword{Biomarker interaction }
\keyword{Likelihood ratio test}
\keyword{Penalized likelihood ratio test}
\keyword{Predictive effect}
\keyword{Prognostic effect}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
