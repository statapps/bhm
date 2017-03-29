\name{pIndexControl}
\alias{pIndexControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary function for pIndex fitting}
\description{
Auxiliary function for \code{\link{pIndex}} fitting. 
Typically only used internally by 'pIndexFit', but may be used to construct a control argument to either function. 
}
\usage{
  pIndexControl(method = c("Efron", "Elc", "Elw"), ci = c("Bootstrap", "Jackknife"), 
                         alpha = 0.05, B = 0, pct = 0.5) 
}

\arguments{
  \item{method}{choose either `Efron' for Efron method, 
   'Elc' for conditional empirical likelihood, 
    or `Elw' for weighted empirical likelihood method. The default value is `Efron'}
  \item{ci}{Method to construct confidence interval, `Bootstrap' for Bootstrap method and 
   `Jackknife' for Jackknife method}
  \item{B}{number of Bootstrap sample}
  \item{alpha}{significance level (e.g. alpha=0.05)}
  \item{pct}{Percentile of threshold (i.e. the cut point), default is 0.5}
}

\details{
Control is used in model fitting of `pIndex'.
}
\value{
This function checks the internal consisitency and returns a list of value as inputed to control model fit of pIndex.
}

\author{Bingshu E. Chen}
\note{
  Based on code from Bingshu E. Chen.
}

\seealso{
\code{\link{bhm}, \link{pIndex}}
}
\examples{
## To calculate the probability index for a biomarker with conditional empirical likelihood method, 
## and the corresponding 90 percent CI using Bootstrap method with 10000 bootstrap sample

ctl = pIndexControl(method = 'Elc', ci = 'Bootstrap', B = 10000, alpha = 0.1)

##
## then fit the following model
##
#  fit = pIndex(y~x1 + x2, family = 'surv', control = ctl)
##
}
\keyword{ control }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
