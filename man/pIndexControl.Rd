\name{pIndexControl}
\alias{pIndexControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary function for pIndex fitting}
\description{
Auxiliary function for \code{\link{pIndex}} fitting. 
Typically only used internally by 'pIndexFit', but may be used to construct a control argument to either function. 
}
\usage{
  pIndexControl(method = c("Efron", "Elc", "Elw", "Pic"), 
                model = c("default", "local", "threshold"), 
	        ci = c("Bootstrap", "Jackknife"), weights = NULL, 
		kernel = NULL, h = 0.1, w  = seq(0.05, 0.95, 0.05), 
		alpha = 0.05, B = 0, pct = 0.5, tau=NULL) 
}

\arguments{
  \item{method}{choose either `Efron' for Efron method, 
   `Elc' for conditional empirical likelihood, 
   `Elw' for weighted empirical likelihood method, and 
   `Pic' for piecewise exponential distribution. The default value is `Efron'}
  \item{model}{`default' for default pIndex model, `local' for kernel method, `threshold' for threshold method}
  \item{ci}{Method to construct confidence interval, `Bootstrap' for Bootstrap method and 
   `Jackknife' for Jackknife method}
  \item{weights}{case weight}
  \item{kernel}{kernel funtion types, including "gaussian", "epanechnikov", "rectangular", "triangular", "biweiht", "cosine", "optcosine". The default value is `gaussian'}

  \item{h}{bandwidth, defaul is 0.1}
  \item{w}{percentile of biomarker value for local fit}
  \item{B}{number of Bootstrap sample}
  \item{alpha}{significance level (e.g. alpha=0.05)}
  \item{pct}{Percentile of threshold (i.e. the cut point), default is 0.5}
  \item{tau}{maximum time tau to be used for pIndex}
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
