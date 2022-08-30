\name{rmscb}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{rmscb}
\alias{rmsControl}
\alias{rmscb.default}
\alias{rmscb.formula}
\alias{plot.rmscb}
\alias{print.rmscb}
\alias{summary.rmscb}
\alias{summary.rmscb}
\alias{K_func}
\alias{getdat}
\alias{getdat1}
\alias{getdat2}

\title{
Fitting Restricted Mean Survival Time Models with a Continuous Biomarker
}
\description{
\{rmscb\} is an R function for restricted mean survival time (RMST) as a continuous function for a 
biomarker variables.
The model is specified by giving a symbolic description of the linear predictor and a description of the distribution family.
}
\usage{
  rmscb(x, ...)

  \method{rmscb}{formula}(formula, data, subset, na.action, tau=5, h=0.2, w0=NULL, 
        sig.level = 0.95, rho = 2,...)
# use 
#          rmscb(y ~ biomarker)             
#
# to fit a prognostic model with biomarker term only 
# 
# use
#          rmscb(y ~ biomarker + trt)
#
# to fit the difference of RMSTs between two treatment groups.
#
# use 
   \method{rmscb}{default}(x, y, control, ...)
#
# to fit a model without formula, where the biomarker shall be in the 
# first dependent variable
}

\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to that class): 
    a symbolic description of the model to be fitted. 
    The details of model specification are given under 'Details'.}
  \item{data}{an optional data frame, list or environment (or object coercible by 'as.data.frame' 
    to a data frame) containing the variables in the model. 
    If not found in data, the variables are taken from environment(formula), 
    typically the environment from which "rmscb" is called.}
  \item{tau}{a prespecified time point at which the restricted mean survival time will be calculated.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. 
    The default is set by the 'na.action' setting of 'options', and is 'na.fail' if that is unset. 
    The 'factory-fresh' default is 'na.omit'.
    Another possible value is 'NULL', no action.  Value 'na.exclude' can be useful.}
  \item{h}{the bandwidth, default is h = 0.2, if h = NULL, then the bandwidth will be selected 
    using the cross validation method.}
  \item{w0}{the values of biomarker at which the RMST E(T|w=w0) will be estimated.}
  \item{sig.level}{the significant level of the simultaneuous confidence band will be constructed.}
  \item{rho}{the mode for the prediction error used in the cross validation bandwidth selection.}
  \item{x}{for 'rmscb.default', x is a design matrix of dimension n * p}
  \item{y}{y is a vector of observations of length n for a "Surv" survival object.}
  \item{control}{a list of parameters for controlling the fitting process. See "rmsControl" for details} 
  \item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below).}  
}

\details{
'biomarker' is a Biomarker variable. This variable is required and shall be the first dependent variable in the formula.

'rmscb.default' is the workhorse functions: they are not normally called directly but can be more efficient
 where the response vector, design matrix and family have already been calculated.
}

\value{
%%  ~Describe the value returned
rmscb returns an object of class inheriting from "rmscb". See later in this section.  

The function "summary" (i.e., "summary.rmscb") can be used to obtain or print a summary of the results, 
for example, the 95 percent CI of the parameters.

An object of class "rmscb" is a list containing at least the following components:
  
  \item{w0}{w0 from the input.}
  \item{rms}{a named vector of restriected mean survival time from the "rmscb".}
  \item{LB}{lower bound of the simultaneuous confidence band.}
  \item{UB}{upper bound of the simultaneuous confidence band.}
}
\references{
Teng, W., Jiang, W. and Chen, B. E. (2022). Continuous threshold models with two-way interactions in survival analysis. 
Statistics in Medicine, 
submitted.
}

\author{
Wen Teng, Wenyu Jiang and Bingshu E. Chen (bingshu.chen@queensu.ca)
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%\code{\link{plot.rmscb}},
%\code{\link{print.rmscb}},
%\code{\link{summary.rmscb}}
%}

\examples{
##
## Generate a random data set
n = 100
age = runif(n, 0, 1)*100
tm = rexp(n, 1/10)
status = rbinom(n, 1, 0.5)
trt= rbinom(n, 1, 0.5)
## fit a restricted mean survival time with one biomarker 

   fit = rmscb(Surv(tm, status)~age)
   print(fit)
## plot(fit)
## summary(fit) 
}

\keyword{Biomarker}
\keyword{Restricted mean survival time}
\keyword{Predictive effect}
\keyword{Prognostic effect}
