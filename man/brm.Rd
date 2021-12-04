\name{brm}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{brm}
\alias{brm.default}
\alias{brm.formula}
%\alias{predict.brm}
\alias{residuals.brm}
\alias{summary.brm}

\title{
Fitting Biomarker Continuous Threshold Models
}
\description{
\{brm\} is an R function for Continuous Threshold Models.
It uses the maximum likehood method (Liu and Chen, 2020) 
to identify a cut-point (thershold parameter) for the biomarker in the Cox proportional hazards model.
The model is specified by giving a symbolic description of the linear predictor and a description of the distribution family.
}
\usage{
brm(x, ...)

\method{brm}{formula}(formula, data = list(...), interaction = TRUE, 
      method = c("gradient", "profile", "single"), q=1, epsilon = NULL, ...)
# use 
#          brm(y ~ biomarker)             
# or
#          brm(y ~ biomarker + x1 + x2, interactin = FALSE)
#
# to fit a prognostic model with biomarker term only (will be implement in the future)
#
# use 
#
#          brm(y ~ biomarker+treatment+x1+x2+...)
#
# to fit a predictive model with interaciton between biomarker 
# and treatment, adjusted for x1, x2, etc.
# 
# use 
#          brm(x, y, control, ...)
#
# to fit a model without formula
#
# Biomarker shall be in the first dependent variable
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to that class): a symbolic description of the model to be fitted. 
  The details of model specification are given under 'Details'.}

  %\item{family}{a description of the response distribution and link function to be used in the model. 
%The available family function are either "binomial" for fitting a logistic regression model or "surv" for fitting a Cox proportional hazards model}

  \item{data}{an optional data frame, list or environment (or object coercible by 'as.data.frame' to a data frame) containing the variables in the model. 
If not found in data, the variables are taken from environment(formula), typically the environment from which glm is called.}

  \item{interaction}{model with interaction term between biomarker and treatment variables. Default is TRUE.}

  \item{method}{method to fit a brm model, such as a gradient method or a single-index model. The default method is the "gradient" method.}

  \item{q}{nubmer of biomarker variables. If q > 1, a single-index model will be fitted. Default is q = 1.}

  \item{x}{For "brm.default", x is a design matrix of dimension n * p and y is a vector of observations of length n for a "Surv" survival object for "coxph" models.}

  %\item{control}{a list of parameters for controlling the fitting process. See "brmControl" for details} 

  \item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below).}  
  \item{epsilon}{Step width for the profile likelihood method, default is (max(w)-min(w))/20.}
}
\details{
'biomarker' is a Biomarker variable. This variable is required and shall be the first dependent variable in the formula.

"interaction" is an option of fitting model with itneractin term.
When interaction = TRUE, a predictive biomarker model will be fitted.
When interaction = FALSE, a prognostic biomarker model will be fitted.
Both Biomarker and Treatment variables are required if 'interaction' = TRUE and 'treatment' shall be the second variable in the formula.

"brm.default" is the workhorse functions: they are not normally called directly but can be more efficient where the response vector, design matrix and family have already been calculated.

}
\value{
%%  ~Describe the value returned
brm returns an object of class inheriting from "brm" which inherits from the class glm or 'coxph'. See later in this section.  

The function "summary" (i.e., "summary.brm") can be used to obtain or print a summary of the results, 
for example, the 95% CI of the parameters.

An object of class "brm" is a list containing at least the following components:
  
  \item{coefficients}{a named vector of coefficients from 'brm'}
  \item{c.max}{the maximum likelihood estimate for the threshold parameter(s).}
  \item{theta}{a vector contains both the coefficients and c.max.}
  \item{var}{the variance matrix of theta.}
  \item{loglik}{the log-likelihood with the final values of the coefficients.}
  \item{linear.predictors}{the vector of linear predictors, one per subject.}
  \item{first}{the first derivative vector at the solution.}
  %\item{se}{estimate of the standard error.}
}
\references{
Liu, S. S. and Chen, B. E. (2020). Continuous threshold models with two-way interactions in survival analysis. Canadian Journal of Statistics. 
}
\author{
Shuoshuo Liu (shuoshuo.liu@psu.edu) and Bingshu E. Chen (bingshu.chen@queensu.ca)
}
%\note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bhm}}, 
\code{\link{coxph}},
\code{\link{plot.brm}},
\code{\link{print.brm}},
\code{\link{residuals.brm}},
\code{\link{summary.brm}},
}
\examples{
##
## Generate a random data set
n = 100
b = c(0.5, 1, 1.5)
data = gendat.surv(n, c0 = 0.40, beta = b, type="brm")
age = runif(n, 0, 1)*100
tm = data[, 1]
status = data[, 2]
trt = data[, 3]
ki67 = data[, 4]
## fit a biomarker threshold survival model with one single cut point

   fit = brm(Surv(tm, status)~ki67+trt+age)
##
## fit a prognostic continuous threshold model with biomarker only
##
#    fit = brm(Surv(tm, status)~ki67)
##
## To print the output, use
##
#    print(fit)
##
}
\keyword{Biomarker interaction }
\keyword{Continuous threshold models }
\keyword{Predictive effect}
\keyword{Prognostic effect}
