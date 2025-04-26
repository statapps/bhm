\name{resboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{resboot}
\alias{resboot.default}
\alias{resboot.formula}

\title{
Rresidual Bootstrap Test (RBT) for treatment-biomarker interaction
}
\description{
\{resboot\} is a function to test the existance of treatment-biomarker 
interaction in biomarker threshold model 

g(Y) = b0+b1*I(w>c) + b2*z + b3*I(w>c)*z.

}

\usage{
resboot(x, ...)

\method{resboot}{formula}(formula, family, data=list(...), B = 100, epsilon = 0.01, ...)
#
###To test the null hypothesis of interaction between treatment variable  
###(define by z) and biomarker variables (define by w) for survival dataa, 
###use:
#
# fit = resboot(Surv(time, status) ~ w + z + w:z)
#
}

\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to 
  that class): a symbolic description of the model to be fitted. 
  The details of model specification are given under 'Details'.}

  \item{family}{default is family = 'Surv' for survival data.}

  \item{data}{an optional data frame, list or environment (or object coercible by 
'as.data.frame' to a data frame) containing the variables in the model. 
If not found in data, the variables are taken from environment(formula), typically the enviro
nment from which resboot is called.}

  \item{x}{Here covariate x is a design matrix of dimension n * 1 (for two sample test) 
or dimension n * 2 (for treatment * biomarker interaction).}

  \item{B}{Number of bootstraps, default is B = 100}

  \item{epsilon}{Biomarker (transformed) step length for profile likelihood method, default is epsilon = 0.01}

  \item{\dots}{additional arguments to be passed to the low level regression fitting 
functions (see below).}  
}

\details{
resboot(y~w + z + w:z) will give residual bootstrap p-value for interaction
between biomarker variable (w) and treatment variable (z).
The null hypothesis is given by H0: b3 = 0, where b3 is the regression 
coefficient for the interaction term I(w>c)*z.
Function print(x) can be used to print a summary of resboot results.
}

\value{
%%  ~Describe the value returned
resboot returns an object of class inheriting from "resboot". 
When B > 0, an object of class "resboot" is a list containing at least the following components:
  \item{theta}{the estimated maximum of likelihood ratio statistics}  
  \item{theta.b}{Bootstrap sample of theta}
  \item{sd}{standard deviation of theta based on resampling}
  \item{ci}{(1-alpha) percent confidence interval for theta based on resampling}
}

\references{
Gavanji, P., Chen, B. E. and Jiang, W.(2018). Residual Bootstrap test for interactions in biomarker threshold models with survival data. Statistics in Biosciences.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}

\note{
    Based on code from Parisa Gavanji.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%\code{\link{prolik}},
\code{\link{bhm}}
\code{\link[survival]{coxph}}
}
\examples{
##
## Generate a random data set
n = 30
b = c(0.5, 1, 1.5)
data = gendat.surv(n, c0 = 0.40, beta = b)
tm = data[, 1]
status = data[, 2]
trt = data[, 3]
ki67 = data[, 4]
#
### No run
# 
# fit = resboot(Surv(tm, status) ~ ki67+trt+ki67:trt) 
#
}
\keyword{ Biomarker interaction }
\keyword{ Predictive effect}
\keyword{ Residual bootstrap }
