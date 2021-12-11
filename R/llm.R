### fit a L-shape linear function
llm = function(x, ...) UseMethod("llm")

llm.default = function(x, y, epsilon, ...) {
  ####possible change points
  #xc = seq(min(x), max(x), (max(x)-min(x))/40)
  x2 = x[, 2]
  xc = quantile(x2, seq(0.025, 0.975, epsilon))
  xc.n = length(xc)
  n = length(x2)
  
  #print(xc)
  rsm = 2*n
  c.max = xc[1]
  for (i in 1:xc.n) {
    ci = xc[i]
    xr = ifelse(x2>ci, (x2-ci), 0)
    X = cbind(x, xr)
    lmf = lm.fit(X, y)
    rsm1 = sum(lmf$residuals^2)
    #print(rsm1)
    if(rsm1 < rsm) {
      fit = lmf
      rsm = rsm1
      c.max = ci
    }
  }
  fit$c.max = c.max
  fit$loglik = -2*rsm
  fit$X = x
  fit$y = y
  class(fit) = c("llm", "lm")
  return(fit)
}

llm.formula = function(formula, data=list(...), epsilon = 0.025,...) {
  
  mf = model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  fit = llm.default(x, y, epsilon)
  fit$terms = mt
  fit$call = match.call()
  return(fit)
}

print.llm = function(x, ...) {
  c0 = x$c.max
  beta = x$coef
  p = length(beta)

  cat("L-shape linear model:\n")
  cat("Cut off value = ", c0, '\n')
  cat("beta = ", beta, '\n')
  cat("Model 1: when x < ", c0, "\n")
  cat("y = ", beta[1])
  if(beta[2]>0) cat ('+') else cat('-')
  cat(abs(beta[2]), '* x')
  
  b20 = beta[1]-beta[p]*c0
  b21 = beta[2]+beta[p]
  
  cat("\nModel 2: when x >= ", c0, "\n")
  cat("y = ", b20)
  if(b21>0) cat ('+') else cat('-')
  cat(abs(b21), '* x')
}

plot.llm = function(x, ...) {
  c0 = x$c.max
  x1 = x$X[, 2]
  x11 = seq(min(x1)-0.5, c0, 0.05)
  x12 = seq(c0, max(x1)+0.5, 0.05)
  y  = x$y
  beta = x$coef
  p = length(beta)
  plot(x1, y)
  #abline(beta[1], beta[2])
  lines(x11, beta[1]+beta[2]*x11)
  abline(v = c0, lty = 2)
  b20 = beta[1]-beta[p]*c0
  b21 = beta[2]+beta[p]
  #abline(b20, b21)
  lines(x12, b20+b21*x12)
}

