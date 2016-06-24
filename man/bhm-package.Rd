\name{bhm-package}
\alias{bhm-package}
\alias{bhm-doc}
\docType{package}
\title{
Biomarker Threshold Models
}
\description{
This package fits a threshold regression model for predictive and prognostic biomarker effect 
for binary data and survival data with an unknown cutoff point of biomarker.
}
\details{
"bhm" is a R package for Biomarker Threshold Models.
Please use the following steps to install 'bhm' package:

1. First, you need to install the 'devtools' package. You can skip this step if you have 'devtools' installed in your R. Invoke R and then type

  install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "bhm" package with R commond

  install_github("statapps/bhm")


"bhm" uses different statistical methods to identify cut-point 
(thershold parameter) for the biomarker in either generalized 
linear models or Cox proportional hazards model.
}
\author{Bingshu E. Chen

Maintainer: Bingshu E. Chen <bingshu.chen@queensu.ca>
}
\references{
Chen, B. E., Jiang, W. and Tu, D. (2014). A hierarchical Bayes model for biomarker subset effects in clinical trials. Computational Statistics and Data Analysis. vol 71, page 324-334.
}
\keyword{bhm}
\seealso{
\code{coxph, }
\code{glm, }
\code{survival}
}
\examples{
# fit = bhm(y~biomarker+treatment)
} 
