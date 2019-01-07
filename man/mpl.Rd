\name{mpl}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{mpl}
\alias{mpl.formula}
\alias{mplFit}

\title{
Joint models for clustered data with robust variance
}
\description{
\{mpl\} is a function to fit Joint model for clustered binary and survival data using maximum penalized likelihood (MPL) with robust Jackknife variance.
}

\usage{
mpl(formula, ...)

\method{mpl}{formula}(formula, formula.cox, formula.cluster, data, weights=NULL,
    subset = NULL, max.iter=100, tol = 0.005, jackknife=TRUE, ...)
#
##use:
#
# fit = mpl(y~x, Surv(time, status) ~ w + z, ~cluster, data = data)
#
}

\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to 
        that class): a symbolic description of the generalized linear model to be fitted for 
        binary data.}

  \item{formula.cox}{an object of class "formula"(or one that can be coerced to
        that class): a symbolic description of the Cox proportiobal hazards model to be fitted
        for survival data.}

  \item{formula.cluster}{an object of class "formula"(or one that can be coerced to
         that class): a symbolic description of the cluster variable.}

  \item{data}{an optional data frame, list or environment (or object coercible by 
'as.data.frame' to a data frame) containing the variables in the model. 
If not found in data, the variables are taken from environment(formula), typically the enviro
nment from which mpl is called.}

  \item{weights}{To be added in the future.}
  \item{subset}{To be added in the future.}

  \item{max.iter}{Maximum number of iterations, default is max.iter = 100}

  \item{tol}{Tolrance for convergence, default is tol = 0.005}

  \item{jackknife}{Jackknife method for variance, default is jackknife = TRUE}

  \item{\dots}{additional arguments to be passed to the low level regression fitting 
functions (see below).}  
}

\details{
mpl(y~x, Surv(time, event)~w+z, ~cluster) will fit penalized likelihood for binary and survival data.
Function print(x) can be used to print a summary of mpl results.
}

\value{
%%  ~Describe the value returned
mpl returns an object of class inheriting from "mpl". 
When jackknife = TRUE, an object of class "mpl" is a list containing the following components:
  \item{theta}{the maximum estimate of the regression coefficients and varaince component}  
  \item{ase}{Asymptotic standard error for theta, which is usually understimated}
  \item{jse}{Jackknife standard error of theta based on resampling, this is considered to be more robust}
}

\references{
Chen, B. E. and Wang, J.(2018). Joint modelling of binary response and survival for clustered data in clinical trials. Statistics in Medicine.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}

\note{
Based on code from J. Wang.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{glm}}
\code{\link{coxph}}
}
\examples{
##
### No run
# 
# fit = mpl(y~x, Surv(tm, status) ~ trt + ki67, ~center.id) 
#
}
\keyword{ Joint Model }
\keyword{ Jackknife }
