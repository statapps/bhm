pIndex = function(x, ...) UseMethod("pIndex") 

pIndex.default = function(x, y, control, ...) {
  cat("pIndex: Porbability Index method for survival data\n")
  x = as.matrix(x)

  xm = as.factor(x[, 1])
  x.n = as.numeric(xm)
  if (max(x.n) > 2) {
     warning("x1 is converted to a binary variable")
     x.n = ifelse(x.n > quantile(x.n, control$pct), 2, 1)
  }
  x[, 1] = x.n - 1
  n = length(y[, 1])
  ci = control$ci
  
  theta = .pIndexFit(x, y, control)
  fit = list(theta = theta)
  B = control$B
  if(ci == "Jackknife") {
    theta.b = rep(0, n)
    for(i in 1:n) {
      x.b = x[-i, ]
      y.b = y[-i, ]
      theta.b[i] = .pIndexFit(x.b, y.b, control)
    }
    fit$theta.b = theta.b
  }
  if(ci == "Bootstrap" & B > 0) {
    theta.b = rep(0, B)
    for(i in 1:B) {
      idx = sample(1:n, n, replace = TRUE)
      x.b = x[idx, ]
      y.b = y[idx, ]
      theta.b[i] = .pIndexFit(x.b, y.b, control)
    }
    fit$theta.b = theta.b
  }
  fit$call = match.call()
  fit$control = control
  if(ncol(x) == 2) fit$interaction = TRUE 
  else fit$interaction = FALSE
  class(fit) = "pIndex"
  return(fit)
}

pIndexControl = function(method = c("Efron", "Elc", "Elw"), ci = c("Bootstrap", "Jackknife"), 
                         alpha = 0.05, B = 0, pct = 0.5) {
  method = match.arg(method)
  ci = match.arg(ci)
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) 
    stop("number of replication 'alpha' must be between 0 and 1")
  if (!is.numeric(B) || B < 0) 
    stop("value of 'B' must be >= 0")
  
  return(list(method = method, ci=ci, alpha = alpha, B = B, pct = pct))
}

pIndex.formula = function (formula, data = list(...), control = list(...), ...) {
  mf = model.frame(formula = formula, data = data)
  x = model.matrix(attr(mf, "terms"), data = mf)
  x = x[, -1]
  y = model.response(mf)
  
  control = do.call("pIndexControl", control)
  fit = pIndex.default(x, y, control, ...)
  fit$call = match.call()
  return(fit)
}

print.pIndex = function(x, ...) {
   control = x$control
   theta = x$theta
   theta.b = x$theta.b
   cat("Call:\n")
   print(x$call)
   sd = sd(theta.b)
   alpha = control$alpha
   ci1 = quantile(theta.b, c(alpha/2, 1-alpha/2))
   if(x$interaction) cat("\nThe probability index for interaction: Pr(T1a<T2a) - Pr(T1b<T2b) \n\n")
   else cat("\nThe probability index for two groups: Pr(T1<T2) \n\n")
   cat("    theta = ", theta, '\n')
   cat("       sd = ", sd, '\n')
   cat("95% confidence interval, \n")
   print(ci1)
}

.pIndexFit = function(x, y, control) {
  x = as.matrix(x)
  x.ncol = ncol(x)
  if(x.ncol>2) stop("x shall be a vector or a matrix with one or two 2 columns.")
  
  if(x.ncol == 2) {
    xm = as.factor(x[, 2])
    x.n = as.numeric(xm)
    if (max(x.n) > 2) {
      warning("x2 is converted to a binary variable")
      x.n = ifelse(x.n > quantile(x.n, control$pct), 2, 1)
    }

    y0 = y[x.n==1, ]
    y1 = y[x.n==2, ]
    x0 = x[x.n==1, 1]
    x1 = x[x.n==2, 1]
  }
  
  method = control$method
  if (method == "Efron") {
    if(x.ncol == 1) theta = .pIndexEfron(x, y)
    if(x.ncol == 2) theta = .pIndexEfron(x1, y1) - .pIndexEfron(x0, y0)
  }
  if (method == "Elc") {
    if(x.ncol == 1) theta = .pIndexElc(x, y)
    if(x.ncol == 2) theta = .pIndexElc(x1, y1) - .pIndexElc(x0, y0)
  }
  if (method == "Elw") {
    if(x.ncol == 1) theta = .pIndexElw(x, y)
    if(x.ncol == 2) theta = .pIndexElw(x1, y1) - .pIndexElw(x0, y0)
  }
  return(theta)
}

### Efron method
.pdfcdf = function(y, group) {
  fit = survfit(y~1)
  cdf = 1-fit$surv
  pdf = diff(c(0, cdf))
  z = rep(group, length(cdf))
  tm = fit$time
  sf = cbind(tm, pdf, z, cdf)
  sf = sf[pdf>0, ]
  return(sf)
}

.boxCox = function(x, r){
  if(r == 0) tx = log(x)
  else tx = (x^r - 1)/r
  return(tx)
}

.pIndexEfron = function(x, y){
  cf1 = .pdfcdf(y[x == 1, ], 1)
  cf2 = .pdfcdf(y[x == 0, ], 0)

  t1 = cf1[, 1]
  t0 = cf2[, 1]
  p1 = cf1[, 2]
  p0 = cf2[, 2]
  
  m1 = length(t1)
  m0 = length(t0)
  
  theta = 0
  #theta1 = 0
  for (i in 1:m0) {
    for (j in 1:m1) {
      dij = ifelse(t0[i]<=t1[j], 1, 0)
      theta = theta + p0[i]*p1[j]*dij
      #theta1 = theta1 + dij
    }
  }  #theta1 = theta1/(n*n)*4
  return(theta)
}

#### Conditional empirical likelihood method ##############
.pc = function(beta, r, m1, m2, W) {
  p = 1/(m1+m2*exp(beta[1]+beta[2]*.boxCox(W, r)))
  return(p)
}

.pc2 = function(beta, r, m1, m2, W) {
  p = exp(beta[1]+beta[2]*.boxCox(W, r))
  return(p)
}

.ellc = function(beta, r, m1, m2, W) {
  m = m1+m2
  p = .pc(beta, r, m1, m2, W)

  W2 = W[(m1+1):m]
  ell = sum(log(p)) + sum(beta[1]+beta[2]*.boxCox(W2, r))
  return(-ell)
}

.pIndexElc = function(x, y) {
  cf1 = .pdfcdf(y[x==1, ], 1)
  cf2 = .pdfcdf(y[x==0, ], 0)
  
  t1 = cf1[, 1]
  t2 = cf2[, 1]
  
  W = c(t1, t2)
  m1 = length(t1)
  m2 = length(t2)
  m = m1 + m2
  
  rx = seq(-2, 3, 0.1)
  lx = rx
  r.max = rx[1]
  beta.max = c(0, 0)
  elm = 1e10
  for (i in 1:length(rx)){
    fitb = try(nlm(.ellc, beta.max, r = rx[i], m1 = m1, m2 = m2, W = W))
    if(class(fitb) == "try-error") next
    ell = fitb$minimum
    if (ell < elm) {
      elm = ell
      r.max = rx[i]
      beta.max = fitb$estimate
    }
    lx[i] = ell
  }

  p = .pc(beta.max, r.max, m1, m2, W)
  p2 = p*.pc2(beta.max, r.max, m1, m2, W)
  theta = 0
  for (i in 1:m) {
    for (j in 1:m){
      dij =  ifelse(W[i]<= W[j], 0, 1)
      theta = theta+p[i]*p2[j]*dij
    }
  }
  return(theta)
}

##### weighted empirical likelihood methdo ##############
.pw = function(beta, r, m1, m2, W) {
  m = m1 + m2
  U = c(rep(m1/m, m1), rep(m2/m, m2))
  u1 = sum(U[1:m1])
  u2 = sum(U[(m1+1):m])
  p = (u1+u2*exp(beta[1]+beta[2]*.boxCox(W, r)))
  return(U/p)
}

.ellw = function(beta, r, m1, m2, W) {
  m = m1+m2
  p = .pw(beta, r, m1, m2, W)
  U = c(rep(m1/m, m1), rep(m2/m, m2))
  
  W2 = W[(m1+1):m]
  U2 = U[(m1+1):m]
  ell = m*sum(U*log(p)) + m*sum(U2*(beta[1]+beta[2]*.boxCox(W2, r)))
  return(-ell)
}

.pIndexElw = function(x, y) {
  cf1 = .pdfcdf(y[x==1, ], 1)
  cf2 = .pdfcdf(y[x==0, ], 0)
  
  t1 = cf1[, 1]
  t2 = cf2[, 1]
  
  W = c(t1, t2)
  m1 = length(t1)
  m2 = length(t2)
  m = m1 + m2
  rx = seq(-2, 3, 0.1)
  lx = rx
  r.max = rx[1]
  beta.max = c(0, 0)
  elm = 1e10
  for (i in 1:length(rx)){
    fitb = try(nlm(.ellw, beta.max, r = rx[i], m1 = m1, m2 = m2, W = W))
    if(class(fitb) == 'try-error') next
    ell = fitb$minimum
    if (ell < elm) {
      elm = ell
      r.max = rx[i]
      beta.max = fitb$estimate
    }
    lx[i] = ell
  }
  
  p = .pw(beta.max, r.max, m1, m2, W)
  p2 = p*.pc2(beta.max, r.max, m1, m2, W)
  
  theta = 0
  for (i in 1:m) {
    for (j in 1:m){
      dij =  ifelse(W[i]<= W[j], 0, 1)
      theta = theta+p[i]*p2[j]*dij
    }
  }
  return(theta)
}
