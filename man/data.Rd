\name{data}
\alias{data}
\alias{glm.gendat}
\alias{surv.gendat}
\docType{data}
\title{dataset}

\description{
dataset for biomarker threshold model (bhm)
}
\usage{surv.gendat(n, c0, beta)
}
\arguments{
  \item{n}{sample size}
  \item{c0}{cut off point, for example c0 = 0.4}
  \item{beta}{regression coefficient, for example, beta = c(0.3, log(0.5), log(0.25))}
}
\format{
  The format of the data set for analysis shall be a data frame with a response variable (either a 
Surv object for Cox model or a glm response variable object) and at least one depend variable as 
the biomarker variable.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
data set of prostate cancer in the 'survival' package is used as an example in paper 
by Chen, et al. (2014).
}
\source{
prosate dataset can be loaded with 'library(survival)'.
}
\references{
Chen, B. E., Jiang, W. and Tu, D. (2014). A hierarchical Bayes model for biomarker subset eff
ects in clinical trials. Computational Statistics and Data Analysis. vol 71, page 324-334.
}
\examples{
#data(data)
## maybe str(data) ; plot(data) ...
c0 = 0.4
b = c(-0.5, 1.5, 1.3)
data = surv.gendat(n=300, c0 = c0, beta = b)
}
\keyword{datasets}
