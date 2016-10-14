\name{pIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{pIndex}
\alias{pIndex.default}
\alias{pIndex.formula}
\alias{pIndexControl}
\alias{pIndexFit}
\alias{print.pIndex}

\title{
Probability Index for Survival Time Difference
}
\description{
\{pIndex\} is a function to estimate and test differce of survival time among groups.
It is defined as p = Pr{ T_1 < T_2 }, where $T_1$ is survival time for subjects in 
group 1 and $T_2$ is surval time in group 2.
}

\usage{
pIndex(x, ...)

\method{pIndex}{formula}(formula, data, control = list(...),...)

pIndexFit(x, y, control)
}

\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to that class): a sy
mbolic description of the model to be fitted. 
  The details of model specification are given under 'Details'.}

  \item{data}{an optional data frame, list or environment (or object coercible by 
'as.data.frame' to a data frame) containing the variables in the model. 
If not found in data, the variables are taken from environment(formula), typically the enviro
nment from which pIndex is called.}

  \item{x, y}{For 'pIndexFit', x is a design matrix of dimension n * 1 (for two sample test) 
or dimension n * 2 (for treatment * biomarker interaction) and 
y is a "Surv" survival object.}

  \item{control}{a list of parameters for controlling the fitting process. 
See 'pIndexControl' for details} 

  \item{\dots}{additional arguments to be passed to the low level regression fitting 
functions (see below).}  
}

\details{
'biomarker' is a Biomarker variable.
print(x) function can be used to print a summary of pIndex results.
}

\value{
%%  ~Describe the value returned
pIndex returns an object of class inheriting from "pIndex". See later in this section.  
}

\references{
Jiang, S., Chen, B. E. and Tu, D.(2016). Inference on treatment-covariate interaction based on a nonparametric measure of treatment effects and censored survival data. Statistics in Medicine. 35, 2715-2725.
}
\author{
Bingshu E. Chen (bingshu.chen@queensu.ca)
}

\note{
The logistic regression part are based on codes wrote by Tian Fang. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% code{\link{prolik}},
\code{\link{bhm}}, 
}
\examples{
##
## Generate a random data set
n = 300
b = c(0.5, 1, 1.5)
data = surv.gendat(n, c0 = 0.40, beta = b)
age = runif(n, 0, 1)*100
tm = data[, 1]
status = data[, 2]
trt = data[, 3]
ki67 = data[, 4]
}
