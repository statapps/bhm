# fit threshold method using profile likelihood method
#prolik = function(x, ...) UseMethod("prolik")

prolikFit = function(x, y, family, control) {
#  cat("prolik: Profile likelihood  model for covariate defined subset effect\n")
  R = control$R
  x = as.matrix(x)
  fit = pro.fit(x, y, family, control)
  cg = NULL
  cqtl = NULL
  if (R > 0) {
    cg = matrix(0, R, control$c.n)
    for (b in 1:R){
      ftb = pro.fit(x, y, family, control)
      cg[b, ] = ftb$c.max
    }
  
    alpha = control$alpha/2
    ptl  = c(alpha, 1-alpha)
    cqtl = apply(cg, 2, quantile, ptl)
  }
  cfit = fit$c.fit
  pfit = list(cg = cg, c.max = fit$c.max, cqtl=cqtl, coefficients = cfit$coefficients, StdErr = sqrt(diag(vcov(cfit))), c.fit = cfit, var_names = colnames(x))
  return(pfit)
}

pro.fit = function(x, y, family, control){
  c.n = control$c.n
  epsilon = control$epsilon

  lik = -nrow(x)*10
  for (u in seq(0.05, 0.95, epsilon)) {
    if (c.n == 2) {
      for (u1 in seq(0.05, 0.9, epsilon)) {
        if (u1 < (u-0.05)) {
          cu = c(u1, u)
          cfit = thm.fit(x, y, family, cu)
          lglk= logLik(cfit)
          if(length(lglk) == 2) lglk = lglk[2]
          if (lglk > lik) {
            lik = lglk
            cx = cu
          }
        }
      }
    } else {
      cfit = thm.fit(x, y, family, u)
      lglk= logLik(cfit)
      if(length(lglk) == 2) lglk = lglk[2]
      if (lglk > lik) {
        lik = lglk
        cx = u
      }
    }
  }
  fit = list(c.max=cx, coefficients=cfit$coefficients, c.fit = cfit)
  return(fit)
}
