pIndex = function(x, ...) UseMethod("pIndex") 

pIndex.default = function(x, y, control, ...) {
  cat("pIndex: Porbability Index method for survival data\n")
  x = as.matrix(x)
  x.ncol = ncol(x)
  if(x.ncol>2) stop("x shall be a vector or a matrix with one or two 2 columns.")
  
  xm = as.factor(x[, 1])
  x.n = as.numeric(xm)
  if (max(x.n) > 2) {
    cat("\nNote: x1 is converted to a binary variable using", control$pct*100, 
        '% percenttile\n', sep = "")
    x.n = ifelse(x.n > quantile(x.n, control$pct), 2, 1)
  }
  x[, 1] = x.n - 1

  # kernel is null, use percentage cutpoint
  if(x.ncol==2 & is.null(control$kernel)) {
    xm = as.factor(x[, 2])
    x.n = as.numeric(xm)
    if (max(x.n) > 2) {
      cat("\nNote: x2 is converted to a binary variable using ", control$pct*100, 
          '% percentile\n', sep = "")
      x.n = ifelse(x.n > quantile(x.n, control$pct), 2, 1)
    }
    x[, 2] = x.n
  }
  
  n = length(y[, 1])
  ci = control$ci
  
  # use local fit when kernel is not null
  if(is.null(control$kernel)) {
    theta = .pIndexFit(x, y, control)
    fit = list(theta = theta)
  } else {
    xc = x
    xc[, 2] = x.cdf(x[, 2])
    fit = .pIndexFitLocal(xc, y, control)
    fit$w0 = quantile(x[, 2], fit$w, type = 3)
  }

  B = control$B
  if(ci == "Jackknife") {
    theta.b = rep(0, n)
    cat("\n\nJackknife resampling...\n")
    if(n > 100){
      cat("|=>                                             Done!\n|")
    } 
    for(i in 1:n) {
      x.b = x[-i, ]
      y.b = y[-i, ]
      theta.b[i] = .pIndexFit(x.b, y.b, control)
      if(n > 100 & i %% floor(n/40) == 0) cat("=")
      control$B = n
    }
    fit$theta.b = theta.b
    theta.bar = mean(theta.b)
    fit$sd = sqrt((n-1)/n*sum((theta.b-theta.bar)^2))
    B = n
  }
  if(ci == "Bootstrap" & B > 0) {
    theta.b = rep(0, B)
    cat("\nBootstrap resampling...\n")
    if(B > 100){
      cat("|=>                                              Done!\n|")
    } 
    for(i in 1:B) {
      idx = sample(1:n, n, replace = TRUE)
      x.b = x[idx, ]
      y.b = y[idx, ]
      theta.b[i] = .pIndexFit(x.b, y.b, control)
      if(B > 100 & i %% floor(B/40) == 0) cat("=")
    }
    fit$theta.b = theta.b
    fit$sd = sd(theta.b)
  }
  cat('\n')
  fit$call = match.call()
  fit$control = control
  z.a = qnorm(1-control$alpha/2)
  if(B>0) fit$ci = c(theta - z.a*fit$sd, theta + z.a*fit$sd)
  if(ncol(x) == 2) fit$interaction = TRUE 
  else fit$interaction = FALSE
  class(fit) = "pIndex"
  return(fit)
}

pIndexControl = function(method = c("Efron", "Elc", "Elw"), ci = c("Bootstrap", "Jackknife"), 
          weights=NULL, kernel = NULL, h=0.1, w=seq(0.05, 0.95, 0.05), alpha = 0.05, B = 0, pct = 0.5) {
  method = match.arg(method)
  ci = match.arg(ci)
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) 
    stop("number of replication 'alpha' must be between 0 and 1")
  if (!is.numeric(B) || B < 0) 
    stop("value of 'B' must be >= 0")

  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")

  if(!is.null(kernel) && B>0) stop("Bootstrap for local estimate coming soon")
  
  return(list(method = method, ci=ci, weights=weights, kernel=kernel, h=h, w=w, alpha = alpha, B = B, pct = pct))
}

pIndex.formula = function (formula, data = list(...), control = list(...), ...) {
  mf = model.frame(formula = formula, data = data)
  x = model.matrix(attr(mf, "terms"), data = mf)
  x = x[, -1]
  y = model.response(mf)
  
  control = do.call("pIndexControl", control)
  fit = pIndex.default(x, y, control, ...)
  fit$call = match.call()
  fit$x = x
  return(fit)
}

print.pIndex = function(x, ...) {
  control = x$control
  theta = x$theta
  theta.b = x$theta.b
  kernel = control$kernel
  cat("Call:\n")
  print(x$call)

  if(x$interaction) cat("\nThe probability index for interaction: Pr(T1a<T2a) - Pr(T1b<T2b) \n\n")
  if(is.null(control$kernel)) cat("\nThe probability index for two groups: Pr(T1<T2) \n\n")

  #print probability index for regular fit
  if(is.null(kernel)) cat("    theta = ", theta, '\n')

  #print list of probability index for local smooth fit
  if(!is.null(kernel)) {
    cat('\nSmooth probability index for Pr(T1<T2|W = w): \n\n')
    out = cbind(w=x$w0, theta = x$theta)
    out = round(out*10000)/10000
    colnames(out) = c(colnames(x$x)[2], 'theta')
    print(out[ceiling(length(x$w0)*seq(0.1, 0.9, 0.1)), ])
  }
  if(control$B > 0) {
    sd = x$sd
    alpha = control$alpha
    ci = quantile(theta.b, c(alpha/2, 1-alpha/2))
    ci[1] = x$ci[1]
    ci[2] = x$ci[2]
    cat("       sd = ", sd, '\n\n', control$ci, '')
    cat((1-control$alpha)*100, "% confidence interval based on B = ", control$B, " resampling\n", sep = "")
    print(ci)
  }
}

.pIndexFit = function(x, y, control) {
  x = as.matrix(x)
  x.ncol = ncol(x)
  if(x.ncol>2) stop("x shall be a vector or a matrix with one or two 2 columns.")

  weights = control$weights
  if(is.null(weights)) 
    weights = rep(1, nrow(x))
  
  if(x.ncol == 2) {
    x.n = x[, 2]
    y0 = y[x.n==1, ]
    y1 = y[x.n==2, ]

    x0 = x[x.n==1, 1]
    x1 = x[x.n==2, 1]

    w0 = weights[x.n==1]
    w1 = weights[x.n==2]
  }
  
  method = control$method
  if (method == "Efron") {
    if(x.ncol == 1) theta = .pIndexEfron(x, y, weights)
    if(x.ncol == 2) theta = .pIndexEfron(x1, y1, w1) - .pIndexEfron(x0, y0, w0)
  }
  if (method == "Elc") {
    if(x.ncol == 1) theta = .pIndexElc(x, y, weights)
    if(x.ncol == 2) theta = .pIndexElc(x1, y1, w1) - .pIndexElc(x0, y0, w0)
  }
  if (method == "Elw") {
    if(x.ncol == 1) theta = .pIndexElw(x, y, weights)
    if(x.ncol == 2) theta = .pIndexElw(x1, y1, w1) - .pIndexElw(x0, y0, w0)
  }
  return(theta)
}

### Kernal fit for local biomarker variable
.pIndexFitLocal=function(x, y, control) {
  h = control$h
  if(is.null(control$w))
    w = sort(x[, 2])
  else w = control$w

  x1 = x[, 1]
  x2 = x[, 2]
  theta = w
  control_w = control
  for(i in 1:length(w)) {
    wg = .K_func(x2, w[i], control$h, control$kernel)
    control_w$weights = wg
    theta[i] = .pIndexFit(x1, y, control_w)
  }
  return(list(theta = theta, w = w))
}

### Efron method
.pdfcdf = function(y, group, weights) {
  fit = survfit(y~1, weights=weights, subset=(weights>0))
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

.pIndexEfron = function(x, y, weights){
  cf1 = .pdfcdf(y[x==1, ], 1, weights[x==1])
  cf2 = .pdfcdf(y[x==0, ], 0, weights[x==0])

  t1 = cf1[, 1]
  t0 = cf2[, 1]
  p1 = cf1[, 2]
  p0 = cf2[, 2]
  
  m1 = length(t1)
  m0 = length(t0)
  
  theta = 0
  #checking using althernative method theta1
  #theta1 = 0
  for (i in 1:m0) {
    for (j in 1:m1) {
      dij = ifelse(t0[i]<=t1[j], 1, 0)
      theta = theta + p0[i]*p1[j]*dij
      #theta1 = theta1 + dij
    }
  } #theta1 = theta1/(n*n)*4
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

.pIndexElc = function(x, y, weights) {
  cf1 = .pdfcdf(y[x==1, ], 1, weights[x==1])
  cf2 = .pdfcdf(y[x==0, ], 0, weights[x==0])
  
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
    fitb = try(optim(beta.max, .ellc, method = "BFGS", r = rx[i], m1 = m1, m2 = m2, W = W), 
               silent = TRUE)
    if(class(fitb) == "try-error") next
    ell = fitb$value
    if (ell < elm) {
      elm = ell
      r.max = rx[i]
      beta.max = fitb$par
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

.pIndexElw = function(x, y, weights) {
  cf1 = .pdfcdf(y[x==1, ], 1, weights[x==1])
  cf2 = .pdfcdf(y[x==0, ], 0, weights[x==0])
  
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
    fitb = try(optim(beta.max, .ellw, method = "BFGS", r = rx[i], m1 = m1, m2 = m2, W = W), 
               silent = TRUE)
    if(class(fitb) == 'try-error') next
    ell = fitb$value
    if (ell < elm) {
      elm = ell
      r.max = rx[i]
      beta.max = fitb$par
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
