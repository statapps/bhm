\name{control}
\alias{bhmControl}
\alias{prolikControl}
\alias{x.cdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary function for bhm fitting}
\description{
Auxiliary function for \code{\link{bhm}} fitting. 
Typically only used internally by 'bhmFit', but may be used to construct a control argument to either function. 
}
\usage{
  bhmControl(method = 'Bayes', interaction, biomarker.main, alpha, 
               B, R, thining, epsilon, c.n, beta0, sigma0)
}

\arguments{
  \item{method}{choose either 'Bayes' for Bayes method with MCMC or 'profile' for profile likelihood method with Bootstrap. The default value is 'Bayes'}
  \item{interaction}{an option of fitting model with itneractin term
     When interaction = TRUE, a predictive biomarker model will be fitted
     When interaction = FALSE, a prognostic biomarker model will be fitted
     The default value is interaction = TRUE.}

  \item{biomarker.main}{include biomarker main effect, default is TRUE}
  \item{B}{set number of burn in}
  \item{R}{set number of replications for Bayes meothd or number of Bootstrap for profile likelihood method}
  \item{thining}{thining parameter for Gibbs samples, default is 2}
  \item{epsilon}{step length for profile likelihood method, default is 0.01}
  \item{alpha}{specify the significance level (e.g. alpha=0.05)}
  \item{c.n}{specify number of threshold (i.e. cut point), default is 1}
  \item{beta0}{initial value for mean of the prior distribution of beta}
  \item{sigma0}{initial value for varaince of the prior distribution of beta}
}

\details{
Conctrol is used in model fitting of bhm.
}
\value{
This function checks the internal consisitency and returns a list of value as inputed to control model fit of bhm.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Bingshu E. Chen and Tian Fang}
\note{
  Based on code from Tian Fang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bhm}}
}
\examples{
## To fit a prognostic model for biomarker with two cut-points, 
## 500 burn-in samples and 10000 Gibbs samples,

ctl = bhmControl(interaction = FALSE, B = 500, R = 10000, c.n = 2)

##
## then fit the following model
##
##  fit = bhmFit(x, y, family = 'surv', control = ctl)
##
}
\keyword{ control }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
