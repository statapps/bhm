pIndex = function(x, ...) UseMethod("pIndex") 

pIndex.default = function(x, y, control, ...) {
  #cat("pIndex: Porbability Index method for survival data\n")
 
  if(is.null(control$cut)) control$cut = as.vector(quantile(y[, 1], seq(0.2, 1.0, 0.2)))
  #.pIndexPic(x1, y, control)
  #
  x = as.matrix(x)
  x.ncol = ncol(x)
  xm = as.factor(x[, 1])
  x.n = as.numeric(xm) - 1
  if (max(x.n) > 2) {
    cat("\nNote: x1 is converted to a binary variable using", control$pct*100, 
        '% percenttile\n', sep = "")
    x.n = ifelse(x.n > quantile(x.n, control$pct), 1, 0)
  }
  x[, 1] = x.n

  # use cutpoint
  if(x.ncol==2) {
    xm = as.factor(x[, 2])
    x.n = as.numeric(xm) - 1
    if (max(x.n) > 2) {
      cat("\nNote: x2 is converted to a binary variable using ", 
	  control$pct*100, '% percentile\n', sep = "")
      x.n = ifelse(x.n > quantile(x.n, control$pct), 1, 0)
    }
    x[, 2] = x.n 
  }
  fit = pIndexFit(x, y, control)

  fit$call = match.call()
  if(ncol(x) == 2) fit$interaction = TRUE 
  else fit$interaction = FALSE
  return(fit)
}

pIndexControl = function(method = c("Efron", "Elc", "Elw", "Pic"), 
			 model = c("default", "local", "threshold"), 
			 ci = c("Bootstrap", "Jackknife"), 
          weights=NULL, kernel = NULL, h=0.1, w=seq(0.05, 0.95, 0.05), alpha = 0.05, B = 0, pct = 0.5) {
  method = match.arg(method)
  model  = match.arg(model)
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
  
  return(list(method = method, model = model, ci=ci, weights=weights, kernel=kernel, h=h, w=w, alpha = alpha, B = B, pct = pct))
}

pIndex.formula = function (formula, data = list(...), control = list(...), ...) {
  mf = model.frame(formula = formula, data = data)
  x = model.matrix(attr(mf, "terms"), data = mf)
  x = x[, -1]
  y = model.response(mf)
  
  x = as.matrix(x)
  x.ncol = ncol(x)
  if(x.ncol>2) stop("x shall be a vector or a matrix with 1 or 2 columns.")

  control = do.call("pIndexControl", control)
  model = control$model

  if(model == "default") fit = pIndex.default(x, y, control)
  else if(model == "local") fit = pIndexLocal(x, y, control)
  else if(model == "threshold") fit = pIndexThreshold(x, y, control)

  fit$call = match.call()
  fit$x = x
  fit$control = control
  class(fit) = "pIndex"
  return(fit)
}

print.pIndex = function(x, ...) {
  control = x$control
  theta = x$theta
  theta.b = x$theta.b
  model = control$model
  alpha = control$alpha
  z.a = qnorm(1-alpha/2)

  cat("Call:\n")
  print(x$call)

  #print probability index for regular fit
  if(model == "default") {
    if(x$interaction) {
      cat("\nThe probability index for interaction: Pr(T1a<T2a) - Pr(T1b<T2b) \n\n")
    } else {
      cat("\nThe probability index for two groups: Pr(T1<T2) \n\n")
    }
    cat("    theta = ", theta, '\n')
  }

  #print list of probability index for local smooth fit
  if(model=="local") {
    cat('\nSmooth probability index for Pr(T1<T2|W = w): \n\n')
    out = cbind(w=x$w0, theta = x$theta)
    out = round(out*10000)/10000
    colnames(out) = c(colnames(x$x)[2], 'theta')
    print(out[ceiling(length(x$w0)*seq(0.1, 0.9, 0.1)), ])
  }

  if((model=="default") & (control$B > 0)) {
    sd = x$sd
    #ci = quantile(theta.b, c(alpha/2, 1-alpha/2))
    ci = c(theta - z.a*sd, theta + z.a*sd)

    cat("       sd = ", sd, '\n\n', control$ci, '')
    cat((1-control$alpha)*100, "% confidence interval based on B = ", control$B, " resampling\n", sep = "")
    print(ci)
  }
  if(model == "threshold") {
    out = cbind(w = x$w, theta = x$theta, sd = x$sd, z = x$theta/x$sd)
    print(out)
  }
}

plot.pIndex = function(x, ...) {
  control = x$control
  model = control$model
  B     = control$B

  z = qnorm(1-control$alpha/2)
  if((model=="threshold")|(model == "local")) {
    w = x$w
    theta = x$theta
    ylim = c(min(theta), max(theta))
    out = cbind(w = x$w, theta = x$theta, sd = x$sd, z = x$theta/x$sd)
    if(B == 0) plot(w, theta, ylim = ylim, type = 'l')
    if(B > 0) {
      sd    = x$sd
      ci1 = theta - z*sd
      ci2 = theta + z*sd
      ylim = c(min(ci1), max(ci2))
      plot(w, theta, ylim = ylim, type = 'l')
      lines(w, ci1, lty = 2)
      lines(w, ci2, lty = 2)
    }
    abline(h = 0)
  }
}
pIndexFit = function(x, y, control) {
  x.ncol = ncol(x)
  B = control$B
  if(is.null(control$weights)) control$weights = rep(1, nrow(x))
  weights = control$weights
 
  theta.sd = NA
  if(x.ncol == 1) {
    theta = .pIndexFit1(x, y, control)
    if(B > 0) theta.sd = .pIndexBoot(x, y, control)$sd
  }
  if(x.ncol == 2) {
    x2 = x[, 2]
    y0 = y[x2==0, ]
    y1 = y[x2==1, ]

    x0 = x[x2==0, 1]
    x1 = x[x2==1, 1]

    w0 = weights[x2==0]
    w1 = weights[x2==1]

    ctl0 = control
    ctl1 = control
    ctl0$weights = w0
    ctl1$weights = w1
    theta = .pIndexFit1(x1, y1, ctl1) - .pIndexFit1(x0, y0, ctl0)
    if(B>0) {
      sd0 = .pIndexBoot(x0, y0, ctl0)$sd
      sd1 = .pIndexBoot(x1, y1, ctl1)$sd
      theta.sd = sqrt(sd0^2 + sd1^2)
    }
  }
  return(list(theta=theta, sd = theta.sd))
}

.pIndexFit1 = function(x, y, control) {
  x = as.matrix(x)
  x.ncol = ncol(x)
  if(x.ncol>1) 
    stop("pIndexFit1: x shall be a vector or a matrix with one columns.")

  weights = control$weights
  method = control$method
  if (method == "Pic")   theta = .pIndexPic(x, y, control)
  if (method == "Efron") theta = .pIndexEfron(x, y, weights)
  if (method == "Elc")   theta = .pIndexElc(x, y, weights)
  if (method == "Elw")   theta = .pIndexElw(x, y, weights)
  return(theta)
}

### Bootstraping, this only works for vector x
.pIndexBoot=function(x, y, control) {
  x = as.matrix(x)
  n = length(y[, 1])
  B = control$B
  ci = control$ci
  fit = NULL

  if(ci == "Jackknife") {
    theta.b = rep(0, n)
    cat("\n\nJackknife resampling...\n")
    if(n > 100){
      cat("|=>                                             Done!\n|")
    }
    for(i in 1:n) {
      x.b = x[-i, ]
      y.b = y[-i, ]
      theta.b[i] = .pIndexFit1(x.b, y.b, control)
      if(n > 100 & i %% floor(n/40) == 0) cat("=")
      control$B = n
    }
    fit$theta.b = theta.b
    theta.bar = mean(theta.b)

    fit$sd = sqrt((n-1)/n*sum((theta.b-theta.bar)^2))
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
      theta.b[i] = .pIndexFit1(x.b, y.b, control)
      if(B > 100 & i %% floor(B/40) == 0) cat("=")
    }
    fit$theta.b = theta.b
    fit$sd = sd(theta.b)
  }
  return(fit)
}

### Kernal fit for local biomarker variable
pIndexLocal=function(x, y, control) {
  xc = x
  xc[, 2] = x.cdf(x[, 2])
  h = control$h
  if(is.null(control$w))
    w = sort(x[, 2])
  else w = control$w
  w0 = quantile(x[, 2], w, type = 3)

  x1 = x[, 1]
  x2 = x[, 2]
  theta = w
  control_w = control
  for(i in 1:length(w)) {
    wg = .K_func(x2, w[i], control$h, control$kernel)
    control_w$weights = wg
    theta[i] = .pIndexFit1(x1, y, control_w)
  }
  return(list(theta = theta, w = w, w0 = w0))
}

pIndexThreshold = function(x, y, control) {
  if(is.null(control$w))
    w = sort(x[, 2])
  else w = control$w

  xc = x
  x2 = x[, 2]
  theta = w
  sd    = w
  control_w = control
  for(i in 1:length(w)) {
    xc[, 2] = ifelse(x2 > w[i], 1, 0)
    fit = pIndexFit(xc, y, control_w)
    theta[i] = fit$theta
    sd[i]    = fit$sd
  }
  fit = list(w = w, theta = theta, sd = sd)
  return(fit)
}

### Efron method
.pdfcdf = function(y, group, weights) {
  fit = survfit(y~1, weights=weights, subset=(weights>0))
  cdf = 1-fit$surv
  pdf = diff(c(0, cdf))
  z  = rep(group, length(cdf))
  tm = fit$time
  sf = cbind(tm, pdf, z, cdf)
  sf = sf[pdf>0, ]
  return(sf)
}

#piecewise constant hazard
.picsf=function(y, group, cut) {
  rownames(y) = c(1:length(y[, 1]))
  dat = data.frame(y)
  spt = survSplit(y~1, data = dat, cut = cut, start = 't0', 
		 end = 't1', episode = 'group')
  spt$group = as.factor(spt$group)

  #total person year
  spt$py = spt$y[, 2] - spt$y[, 1]
  spt$event = spt$y[, 3]
  tab = aggregate(spt[c("event", "py")], by = list(group = spt$group), sum)
  tab$rate = tab$event/tab$py
  return(tab)
}

#mix

pictail=function(y, group, cut) {
  weights = 1
  sf = .pdfcdf(y, group, weights)
  
}

# pIndex for piecewise
.pIndexPic=function(x, y, control){
  cut = control$cut
  #if (is.null(cut)) cut = c(1.5, 2.3, 3, max(y[, 1]))
  sf1 = .picsf(y[x==1, ], 1, cut) 
  sf2 = .picsf(y[x==0, ], 0, cut)
  r1 = sf1$rate
  r2 = sf2$rate
  event = sf1$event + sf2$event
  event = event/(sum(event))

  K = length(cut)
  # theta = sum(r2/(r1+r2)*(1 - exp(-(r1+r2)*(b-a)))*s1(a)*s2(a)))
  rt = r1 + r2
  rx = r2/rt

  #cut = c(cut, cut[K]*10)

  tau = c(0, cut)
  dtau = diff(tau)
  dtau1 = c(0, dtau)
  
  dtau1= dtau1[1:K]

  s0 = exp(-cumsum(rt*dtau1))
  wk = (1 - exp(-rt*dtau))*s0
  wk = event

  theta = sum(rx*wk)
  return(theta)
}

.pIndexEfron = function(x, y, weights){
  cf1 = .pdfcdf(y[x==1, ], 1, weights[x==1])
  cf2 = .pdfcdf(y[x==0, ], 0, weights[x==0])
  t1 = cf1[, 1]
  t0 = cf2[, 1]
  p1 = cf1[, 2]
  p0 = cf2[, 2]
  F0 = cf2[, 4]

  m1 = length(t1)
  m0 = length(t0)
  theta = 0
  #checking using althernative method theta1
  #theta1 = 0
  #for (i in 1:m0) {
  #  for (j in 1:m1) {
  #    dij = ifelse(t0[i]<=t1[j], 1, 0)
  #    theta = theta + p0[i]*p1[j]*dij
  #  }
  #}
  for(i in 1:m1) {
    idx = sum(ifelse(t0 < t1[i], 1, 0))
    if(idx>0) theta = theta + F0[idx]*p1[i]
  }
  #cat('theta = ', theta, theta1, '\n')
  return(theta)
}

#### Conditional empirical likelihood method ##############
.boxCox = function(x, r){
	  if(r == 0) tx = log(x)
  else tx = (x^r - 1)/r
    return(tx)
}

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
