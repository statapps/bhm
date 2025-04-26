\name{bhm-package}
\alias{bhm-package}
\alias{bhm-doc}
\docType{package}
\title{
An R package for the biomarker threshold models
}
\description{
This package fits biomarker threshold regression models for predictive and prognostic biomarker effects 
with binary data and survival data with an unknown biomarker cutoff point (Chen et al, 2014)<DOI:10.1016/j.csda.2013.05.015>. 
Multivariable models can also be fitted for adjusted biomarker effect (Fang et al, 2017)<DOI:10.1016/j.csda.2017.02.011>. 
Tools such as Probability index are included to measure treatment effect, biomarker effect or treatment-biomarker interaction(Jiang et al, 2016)<DOI:10.1002/sim.6907>. 
}
\details{
"bhm" is a R package for Biomarker Threshold Models.
Please use the following steps to install the most recent version of 'bhm' package:

1. First, you need to install the 'devtools' package. 
   You can skip this step if you have 'devtools' installed in your R. Invoke R and then type

   install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "bhm" package from github with R commond

  install_github("statapps/bhm")


"bhm" uses different statistical methods to identify cut-point (thershold parameter) for the 
biomarker in either generalized linear models or Cox proportional hazards model.

A stable version of View the "bhm" package is also available from the Comprehensive R Archive Network 
(https://CRAN.R-project.org/package=bhm) and can be installed using R command 

  install.packages("bhm")

}
\author{Bingshu E. Chen, Tian Fang, Jia Wang, Shuoshuo Liu

Maintainer: Bingshu E. Chen <bingshu.chen@queensu.ca>
}
\references{
Chen, B. E., Jiang, W. and Tu, D. (2014). A hierarchical Bayes model for biomarker subset effects in clinical trials. Computational Statistics and Data Analysis. vol 71, page 324-334.

Fang, T., Mackillop, W., Jiang, W., Hildesheim, A., Wacholder, S. and Chen, B. E. (2017). A Bayesian method for risk window estimatin with application to HPV vaccine trial. Computational Statistics and Data Analysis. 112, page 53-62.

Jiang, S., Chen, B. E. and Tu, D.(2016). Inference on treatment-covariate interaction based on a nonparametric measure of treatment effects and censored survival data. Statistics in Medicine. 35, 2715-2725.

Gavanji, P., Chen, B. E. and Jiang, W.(2018). Residual Bootstrap test for interactions in biomarker threshold models with survival data. Statistics in Biosciences.

Chen, B. E. and Wang, J.(2020). Joint modelling of binary response and survival for clustered data in clinical trials. Statistics in Medicine. Vol 39. 326-339.

Liu, S. S. and Chen, B. E. (2020). Continuous threshold models with two-way interactions in sur
vival analysis. Canadian Journal of Statistics. 

}
\keyword{Biomarker threshold models}
\keyword{Continuous threshold models}
\keyword{Probability index}
\keyword{Indicator threshold models}
\keyword{Residual Bootrap}

\seealso{
\code{\link{bhm}},  
\code{\link{brm}}, 
\code{\link[survival]{coxph}}, 
\code{\link{glm}} 
}
\examples{
# fit = bhm(y~biomarker+treatment)
# print(summary(fit))
} 
