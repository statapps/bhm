\name{bhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{bhm}
\alias{bhm.default}
\alias{bhm.formula}
\alias{bhmFit}
\alias{bhmGibbs}
\alias{prolikFit}
\alias{summary.bhm}
\alias{thm.fit}
\alias{thm.lik}
\alias{x.cdf}

\title{
Fitting Biomarker Threshold Models
}
\description{
\{bhm\} is a R package for Biomarker Threshold Models.
It uses either Hierarchical Bayes method or proflie likehood method (Chen, et al, 2014 and Tian, et al, 2017) 
to identify a cut-point (thershold parameter) for the biomarker in either generalized 
linear models or Cox proportional hazards model.
The model is specified by giving a symbolic description of the linear predictor and a description of the distribution family.
}
\usage{
bhm(x, ...)

\method{bhm}{formula}(formula, family, data, control = list(...),...)

# use 
#          bhm(y ~ biomarker)             
#
# to fit a prognostic model with biomarker term only
#
# use 
#
#          bhm(y ~ biomarker+treatment)
#
# to fit a predictive model with interaciton between biomarker 
# and treatment, use 
#
          bhmFit(x, y, family, control)  
#
# to fit a model without the formula
#
# Biomarker shall be in the first dependent variable
# 
#  To summary a "bhm" boject,
#
          \method{summary}{bhm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to that class): a symbolic description of the model to be fitted. 
  The details of model specification are given under 'Details'.}

  \item{family}{a description of the response distribution and link function to be used in the model. 
The available family function are either "binomial" for fitting a logistic regression model or "surv" for fitting a Cox proportional hazards model}

  \item{data}{an optional data frame, list or environment (or object coercible by 'as.data.frame' to a data frame) containing the variables in the model. 
If not found in data, the variables are taken from environment(formula), typically the environment from which glm is called.}

  \item{x, y}{For "bhmFit", x is a design matrix of dimension n * p and y is a vector of observations of length n for "glm" models or a "Surv" survival object for "coxph" models.}

  \item{control}{a list of parameters for controlling the fitting process. See "bhmControl" for details} 

  \item{object}{object returned from model fit}

  \item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below)}  
}
\details{
'biomarker' is a Biomarker variable. This variable is required and shall be the first dependent variable in the formula.

"interaction" is an option of fitting model with itneractin term.
When interaction = TRUE, a predictive biomarker model will be fitted.
When interaction = FALSE, a prognostic biomarker model will be fitted.
Both Biomarker and Treatment variables are required if 'interaction' = TRUE and 'treatment' shall be the second variable in the formula.

"bhmFit" and "bhmGibbs" are the workhorse functions: they are not normally called directly but can be more efficient where the response vector, design matrix and family have already been calculated.

"x.cdf" is a function that maps biomarker values to interval (0, 1) using its empirical cumulative distribution function. After the threshold parameters are identified, the biomarker variable will be transformed back to its original scale.
}
\value{
%%  ~Describe the value returned
bhm returns an object of class inheriting from "bhm" which inherits from the class glm or 'coxph'. See later in this section.  

The function "summary" (i.e., "summary.bhm") can be used to obtain or print a summary of the results, 
for example, the 95% CI of the parameters.

An object of class "bhm" is a list containing at least the following components:
  
  \item{c.max }{a vector of the mean estimates for the threshold parameter(s)}
  \item{coefficients}{a named vector of coefficients from 'bhm'}
  \item{c.fit}{fitted conditional regression model given c = c.max}
  \item{cg}{Gibbs sample for threshold parmeter c}
  \item{bg}{Gibbs sample for the coefficients beta}
}
\references{
Chen, B. E., Jiang, W. and Tu, D. (2014). A hierarchical Bayes model for biomarker subset effects in clinical trials. Computational Statistics and Data Analysis. vol 71, page 324-334.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}
\note{
The logistic regression part are based on codes wrote by Tian Fang. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% code{\link{prolik}},
\code{\link{brm}}, 
\code{\link[survival]{coxph}}, 
\code{\link{glm}}, 
\code{\link{glmpLRT}}, 
\code{\link{mpl}}, 
\code{\link{pIndex}}, 
\code{\link{resboot}}, 
\code{\link{rmscb}}, 
\code{\link{bhmControl}}
}
\examples{
##
## Generate a random data set
n = 300
b = c(0.5, 1, 1.5)
data = gendat.surv(n, c0 = 0.40, beta = b)
age = runif(n, 0, 1)*100
tm = data[, 1]
status = data[, 2]
trt = data[, 3]
ki67 = data[, 4]
## fit a biomarker threshold survival model with one single cut point

# fit = bhm(Surv(tm, status)~ki67+trt+age, interaction = TRUE, B=5, R=10)


## here B=5 and R=10 is used for test run. In general, B > 500 and R > 2000 is
## recommend for the analysis of biomarker variable. To fit a model with 
## two cut points, use: 
##
##     fit = bhm(Surv(tm, status)~bmk+trt+age, B = 500,  R = 2000, c.n = 2)
##
## To print the output, use
##
#       print(fit)
##
}
\keyword{Biomarker interaction}
\keyword{Indicator threshold models}
\keyword{Logistic regression}
\keyword{Predictive effect}
\keyword{Prognostic effect}
