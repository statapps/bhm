\name{coxpLRT}
\alias{coxpLRT}
\alias{coxpLRT.formula}
\alias{coxpLRT.default}
\alias{resbootCox}
\alias{coxpLRTcvPredErr}

\title{Penalized Likelihood Ratio Test for a Biomarker Threshold Effect in the Cox Proportional Hazards Model}
\description{Fits a Cox proportional hazards model with a smoothed threshold
(change-point) effect for one or more biomarkers, optionally
allowing the threshold effect to have interactions effects with other covariates.
Two methods are developed to test for the presence of the threshold
effect: a penalized likelihood ratio test (LRT) and a residual
bootstrap test. The penalized LRT is based on the asymptotic
distribution of the LRT statistic and provides a fast
approach for hypothesis testing. The residual bootstrap method
approximates the null distribution of the test statistic through
resampling and generally requires substantially more computation time
to run, but provides an alternative when a resampling-based inference
procedure is preferred.

The model is specified by giving a symbolic description of the linear predictor and 
a description of the distribution family.}
\usage{
coxpLRT(XZ, \dots)

\method{coxpLRT}{formula}(formula, formula2 = NULL, formulaList = NULL,
             data = list(...), c0 = 0.5, lambda = 25, K = 50,
             x = FALSE, bootstrap = 0, \dots)

\method{coxpLRT}{default}(XZ, y, lambda = 25, naive.test = FALSE, \dots)
# use 
#          coxpLRT(y ~ biomarker, ~z1+z2, formulaList = c(~1))
# or
#
# to fit a prognostic model with biomarker term adjust for 
# covariates z1, z2, etc.
#
# use 
#
#          coxpLRT(y ~ x1+x2, ~z1+z2+z3+z4+z5, c(~z1, ~z2), ...)
#
# to fit a predictive model with interaciton between biomarker x1 and covariate z1, 
# biomarker x2 and covariate z2, adjusted for z3, z4, z5, etc.
# 
# use 
#          coxpLRT(ZX, y, control, ...)
#
# to fit a model without formula, where ZX = getZX(formula, formula2, formulaList)
#
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the biomarker model to be fitted. The formula cannot be missing. The details of model specification are given under \sQuote{Details}.}
  \item{formula2}{an object of class \code{"formula"}: a symbolic description of the covariate (other than the biomarkers) model to be fitted. The default is \code{NULL}, which means only the biomarker prognostic model will be fitted. The details of model specification are given under \sQuote{Details}.}
  \item{formulaList}{an object of class \code{"formula"}: a symbolic description of the covariate(s) that interact with the biomarkers. The details of model specification are given under \sQuote{Details}.}
  \item{data}{an optional data frame, list, or environment containing the variables in the model.}
  \item{c0}{the cut point for the threshold function. Default is \code{c0 = 0.5}. Other options include \code{c0 = 0.25} and \code{c0 = 0.75}.}
  \item{lambda}{the tuning parameter for the penalized likelihood function.}
  \item{K}{smoothing parameter for the indicator function; the larger the value of \code{K}, the better the approximation of the indicator function. Default is \code{K = 50}.}
  \item{x}{\code{x = TRUE} to return the \code{X} matrix.}
  \item{naive.test}{Conduct a naive LRT test without the penalty term.}
  \item{bootstrap}{
    Integer number of residual bootstrap replicates used to
    approximate the null distribution of the likelihood ratio
    statistic and obtain a bootstrap p-value. If \code{0} (the
    default), no bootstrap is performed and \code{fit$boot$pValue} is
    set to \code{1}.
  }
  \item{XZ}{
    For the default method: the pre-built design object combining the
    biomarker terms (\code{formula}), the other covariates
    (\code{formula2}), and the biomarker-interacting covariates
    (\code{formulaList}), as constructed internally by \code{.getXZ}.
    Not typically supplied directly by the user.
  }
  \item{y}{For the default method: a \code{\link[survival]{Surv}} object of survival times and event/censoring indicators, sorted in decreasing order of time.}
  \item{\dots}{additional arguments to be passed to the low-level regression fitting functions (see below).}
}
\details{
  \code{coxpLRT} is generic, with a formula method
  (\code{coxpLRT.formula}) intended for direct use, and a default
  method (\code{coxpLRT.default}) that operates on an already-built
  design object and performs the actual model fitting and likelihood
  ratio calculation. "coxpLRT.default" is the workhorse functions: they are not normally
called directly but can be more efficient where the response 
vector and design matrix have already been calculated.


  Three formulas jointly specify the model:
  \itemize{
    \item \code{formula} specifies the survival outcome (a
      \code{\link[survival]{Surv}} object on the left-hand side) and
      the biomarker(s) whose threshold (change-point) effect is being
      modeled, e.g. \code{Surv(time, status) ~ x1 + x2}.
    \item \code{formula2} optionally specifies additional covariates
      to be included in the model linearly (i.e. without a threshold
      effect), e.g. \code{~ z}. If omitted, only the biomarker
      prognostic model is fitted.
    \item \code{formulaList} optionally specifies covariates that
      interact with (modify) the biomarker threshold effect(s), e.g.
      \code{c(~z1, ~z2)} for a list of interaction formulas
      corresponding to the biomarkers in \code{formula}.
  }

  Internally, the formula method:
  \enumerate{
    \item Builds the model frame from \code{formula} and extracts the
      survival response, which must be a \code{Surv} object.
    \item Sorts the data in decreasing order of survival/censoring
      time (required by the internal partial likelihood
      calculations).
    \item Constructs the design object via the internal helper
      \code{.getXZ}, combining the biomarker terms, the linear
      covariates from \code{formula2}, and the interaction structure
      from \code{formulaList}, together with the starting cut point
      \code{c0} and the indicator-smoothing parameter \code{K}.
    \item Calls \code{coxpLRT.default} to fit the (penalized) model,
      using \code{lambda} as the penalty tuning parameter, and to
      compute the likelihood ratio statistic for the threshold
      effect.
    \item If \code{bootstrap > 0}, calls the internal helper
      \code{.resbootCox} to generate \code{bootstrap} residual
      bootstrap replicates and attaches the result (including the
      bootstrap p-value) as \code{fit$boot}. Otherwise
      \code{fit$boot$pValue} is set to \code{1}.
  }

  The threshold (cut point) location is a nuisance parameter that is
  unidentified under the null hypothesis of no threshold effect (the
  Davies problem), which is why a bootstrap p-value (via
  \code{bootstrap > 0}) is generally preferred over any asymptotic
  reference distribution for the likelihood ratio statistic returned
  here.
  
  "print", "plot" and "summary" methods can be applied to a fitted "coxpLRT" class to display the results.
}
\value{
  A fitted \code{coxpLRT} object (list), with components returned by
  \code{coxpLRT.default} plus, if requested, an \code{x} component
  (the design matrix) and a \code{boot} component. Typical components
  include:
  \item{coefficients}{Estimated regression coefficients, including the biomarker threshold effect(s) and any interaction terms.}
  \item{cutpoint}{Estimated threshold (cut point) location(s) for the biomarker(s).}
  \item{loglik}{(Penalized) log partial likelihood(s) under the null and/or fitted model, used to form the likelihood ratio statistic.}
  \item{statistic}{The likelihood ratio test statistic for the threshold effect.}
  \item{x}{The model design (\code{X}) matrix, if \code{x = TRUE}.}
  \item{boot}{
    A list with the results of the residual bootstrap, including
    \code{pValue} (the bootstrap p-value; \code{pValue = NA} if
    \code{bootstrap = 0}) and the vector of bootstrap replicate
    statistics for the LRT.
  }
}
\note{
  \code{.getXZ} and \code{.resbootCox} are internal, unexported
  helper functions and are not intended to be called directly by
  users.
}
\references{
  Davies, R. B. (1987). Hypothesis testing when a nuisance parameter
  is present only under the alternative. \emph{Biometrika}, 74(1), 33--43.

  Jiang, W. and Chen, B. E. (2026). Tests for covariate threshold effects for survival data based on penalized Cox regression model with application to biomarker study. Biometrics. Revised manuscripted submitted. 
}
\seealso{
  \code{\link[bhm]{glmpLRT}}, 
  \code{\link{glm}},
  \code{\link[bhm]{plot.glmpLRT}},
  \code{\link[bhm]{print.glmpLRT}}
  \code{\link[survival]{coxph}}, 
  \code{\link[survival]{Surv}}
}
\examples{
  \dontrun{
  ## A simulated example
  bmk = rnorm(100, 3, 0.25)
  age = rnorm(100, 50, 20)
  trt = rbinom(100, 1, 0.5)
  y = Surv(rexp(100), rbinom(100, 1, 0.8))
  dat = data.frame(y, trt, age, bmk)
  #fit = coxpLRT(y~bmk, ~trt+age, c(~bmk), data = dat)
  #print(fit)

## biomarker threshold model, no additional covariates, no interaction
# fit <- coxpLRT(Surv(time, status) ~ x1 + x2,
#               data = dat, c0 = 0.5, lambda = 25, K = 50, bootstrap = 20)
# fit$statistic
# fit$boot$pValue

## with additional (non-biomarker) covariates and biomarker interactions
# fit2 <- coxpLRT(Surv(time, status) ~ x1 + x2, ~ z, c(~z, ~z),
#                data = dat, c0 = c0, lambda = 40, K = K, bootstrap = B)
 }
}
\keyword{Biomarker interaction }
\keyword{Likelihood ratio test}
\keyword{Penalized likelihood ratio test}
\keyword{Predictive effect}
\keyword{Prognostic effect}
\keyword{survival}
%\keyword{htest}
