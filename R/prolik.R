# fit threshold method using profile likelihood method

prolikFit = function(x, y, family, control) {
  R = control$R
  x = as.matrix(x)
  var_names = colnames(x)
  if (family == "surv") var_names = var_names[-1]
  n = length(x[, 1])
  fit = .profit(x, y, family, control)
  cg = NULL
  cqtl = NULL
  if (R > 0) {
    cg = matrix(0, R, control$c.n)
    for (b in 1:R){
      idx = sample(1:n, n, replace = TRUE)
      x.b = x[idx, ]
      if(is.vector(y)) y.b = y[idx]  #For binomial or continuous
      else y.b = y[idx, ] #For survival
      ftb = .profit(x.b, y.b, family, control)
      cg[b, ] = ftb$c.max
    }
  
    alpha = control$alpha/2
    ptl  = c(alpha, 1-alpha)
    cqtl = apply(cg, 2, quantile, ptl)
  }
  cfit = fit$c.fit
  pfit = list(cg = cg, c.max = fit$c.max, cqtl=cqtl, coefficients = cfit$coefficients, StdErr = sqrt(diag(vcov(cfit))), c.fit = cfit, var_names = var_names)
  return(pfit)
}

.profit = function(x, y, family, control){
  c.n = control$c.n
  epsilon = control$epsilon

  lik = -nrow(x)*10
  lglk = lik
  for (u in seq(0.05, 0.95, epsilon)) {
    if (c.n == 2) {
      for (u1 in seq(0.05, 0.9, epsilon)) {
        if (u1 < (u-0.05)) {
          cu = c(u1, u)
          cfit = thm.fit(x, y, family, cu)
          if(cfit$converged) lglk= logLik(cfit)
          if(length(lglk) == 2) lglk = lglk[2]
          if (lglk > lik) {
            lik = lglk
            cx = cu
          }
        }
      }
    } else {
      cfit = thm.fit(x, y, family, u)
      if(cfit$converged) lglk= logLik(cfit)
      if(length(lglk) == 2) lglk = lglk[2]
      if (lglk > lik) {
        lik = lglk
        cx = u
      }
    }
  }
  cfit = thm.fit(x, y, family, cx)
  fit = list(c.max=cx, coefficients=cfit$coefficients, c.fit = cfit)
  return(fit)
}
