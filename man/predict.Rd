\name{predict}
\alias{predict.coxpLRT}
%\alias{predict.glmpLRT}
\title{Predict Method for Model Fits in the \code{bhm} package}
\description{
  Obtains linear predictors from a fitted 
  object, either the linear predictors from the original model fit or
  predictions for new data.
}
\usage{
\method{predict}{coxpLRT}(object, newdata = NULL, \dots)
}
\arguments{
  \item{object}{a fitted model object.}
  \item{newdata}{an optional data frame containing new data for which to compute predictions. If omitted (the default), the linear predictors from the original model fit are returned.}
  \item{\dots}{additional arguments, currently unused, for compatibility with the generic \code{\link[stats]{predict}} method.}
}
\details{
  If \code{newdata} is \code{NULL}, the fitted linear predictors
  stored in \code{object$linear.predictors} are returned directly.
  For most threshold model, these were computed at the estimated threshold location
  \code{object$c.max}. 

  If \code{newdata} is supplied, predictions are instead obtained by
  calling \code{\link[stats]{predict}} on the fitted model object, 
  applied to \code{newdata}.
}
\value{
  A numeric vector of linear predictors, either taken from the
  original \code{coxpLRT} fit (\code{newdata = NULL}) or computed for
  \code{newdata} via the underlying Cox model fit.
}
\note{
  Because the underlying model was fitted using the biomarker
  threshold structure estimated at \code{object$c.max} (and any
  covariates or interaction terms, 
  \code{newdata} must contain columns matching that design
  \emph{exactly} as expected by the stored mode. 
}
\seealso{
  \code{\link{coxpLRT}}, \code{\link[stats]{predict}}
}
\examples{
\dontrun{
#fit <- coxpLRT(Surv(time, status) ~ x1 + x2,
#               data = dat, c0 = 0.5, lambda = 25, K = 50)

## linear predictors from the original fit
#lp <- predict(fit)

## linear predictors for new data
#lp.new <- predict(fit, newdata = newdat)
}
}
\keyword{survival}
\keyword{methods}
