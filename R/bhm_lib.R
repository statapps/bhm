#generate regression coefficients
thm.fit <-function(x, y, family, cx){
  #Remove the intercept term
  c.n = length(cx)
  n.col = length(x[1, ])
  x = x[, -1]
  #remove intercept x1 biomarker, x2 trt 
  if(n.col == 2) x = matrix(x, ncol = 1)
  w = x[, 1]
  if (c.n == 1) x[, 1] = ifelse(cx <= w, 1, 0)
  if (c.n == 2) x[, 1] = ifelse((cx[1]<=w) & (w<=cx[2]), 1, 0)
  
  if (length(grep(":", colnames(x)[3])) == 1) 
    x[, 3] = x[, 1]*x[, 2]
  else if(length(grep(":", colnames(x)[1])) == 1) 
    x[, 1] = x[, 1]*x[, 2]

  x.c_ = x
  if (family == "surv") {
    fit = tryCatch(coxph(y~x.c_), warning = function(w) w)
    if(is(fit, "warning")) fit = list(converged = FALSE)
    else fit$converged = TRUE
  }
  else fit = glm(y~x.c_, family=family)
  
  return(fit)
}

#calculate the joint probability (log likelihood) of y
thm.lik = function(x, y, family, beta, q, cx, control){
  c.n = control$c.n
  beta0 = control$beta0
  s0.inv = control$sigma0.inv
  n.col = length(x[1, ])

  w = x[, 2]
  if (c.n == 1) x[, 2] = ifelse(cx <= w, 1, 0)
  if (c.n == 2) x[, 2] = ifelse((cx[1]<=w) & (w<=cx[2]), 1, 0)

  # x1 intercept, x2 biomarker, x3 trt
  if(length(grep(":", colnames(x)[4])) == 1) 
    x[, 4] = x[, 2]*x[, 3]
  else if(length(grep(":", colnames(x)[2])) == 1)  #no bmk main effect, use x2 for 
    x[, 2] = x[, 2]*x[, 3]                         #interaction

  if (family == "surv") {
    x = x[, -1]
    if(n.col == 2) x = as.matrix(x, ncol = 1)
    xbeta = x%*%beta
    exb = exp(xbeta)
    cxb = log(cumsum(exb))
    lik = sum(y[, 2]*(xbeta-cxb))
  } 
    
  if (family == "binomial") {
    xbeta = x%*%beta
    exb = exp(xbeta)  
    lik = sum(y*xbeta - log(1+exb))
  }

  if (family == "gaussian") {
    xbeta = x%*%beta
    lik = sum((y-xbeta)^2)
  }
  
  if (c.n==1) 
    lik2=dbeta(cx, 2, q, log = TRUE)
  else
    lik2=log(cx[1])+(q[1]-1)*log(cx[2]-cx[1])-q[1]*log(cx[2])+(q[2]-1)*log(1-cx[2])

  if (length(beta) == 1) phib0 = 0.5*(beta-beta0)^2*s0.inv
  else phib0 = 0.5*t(beta-beta0)%*%s0.inv%*%(beta-beta0)

  lik = lik + lik2 - phib0
  return(lik)
}

surv.gendat =  function(n, c0, beta){
  n1 = n/2
  #x = runif(n, -2, 2)
  x = rnorm(n, 0.2, 1)
  x1 = ifelse(x > c0, 1, 0)
  z = c(rep(0, n1), rep(1, n1))
  zx = z*x1
  x0 = rep(1, n)
  X = cbind(z, x1, zx)
  h0 = 0.5
  h = h0*exp(X%*%beta)
  stime = rexp(n, h)               #Failure time.
  endstudy = runif(n, 0, 5)
  cstatus = ifelse(stime>endstudy, 0, 1) ##Censored at end of study time.
  cat('\nCensoring: ', 1-mean(cstatus), '\n')
  stime = ifelse(stime>endstudy, endstudy, stime)
  dat = cbind(time=stime, status=cstatus, z=z, x=x)
  return(dat)
}

glm.gendat =  function(n, c0, beta){
  n1 = n/2
  x = runif(n, 0, 1)
  c.n = length(c0)
  
  if (c.n == 1) x1 = ifelse(x <= c0, 1, 0)
  if (c.n == 2) x1 = ifelse((c0[1] <= x) & (x <= c0[2]), 1, 0)

  z = c(rep(0, n1), rep(1, n1))
  zx = z*x1
  x0 = rep(1, n)
  X = cbind(x0, z, x1, zx)

  eb = exp(X%*%beta)
  p = eb/(1+eb)
  y = rbinom(n, 1, p)
  dat = cbind(y, z, x)
  return(dat)
}

#empirical cdf tranformation used to convert time interval
x.cdf = function(x){
  n = length(x)
  p = rep(0, n)
  for (i in 1:n) {
    p[i] = sum(x<=x[i])
  }
  p = (p-0.5)/n
  return(p)
}

## Kernel function
.K_func = function(w, u, h, kernel = c("gaussian", "epanechnikov", "rectangular", 
		     "triangular", "biweight", "cosine", "optcosine")) {
  kernel = match.arg(kernel)
  x = w-u
  ax = abs(x)
  esp = 1e-40

  kh = switch(kernel, gaussian = dnorm(x, sd = h/2),
    rectangular = ifelse(ax < h, 0.5/h, esp), 
    triangular = ifelse(ax < h, (1 - ax/h)/h, esp),
    epanechnikov = ifelse(ax < h, 3/4 * (1 - (ax/h)^2)/h, esp),
    biweight = ifelse(ax < h, 15/16 * (1 - (ax/h)^2)^2/h, esp),
    cosine = ifelse(ax < h, (1 + cos(pi * x/h))/(2*h), esp),
    optcosine = ifelse(ax < h, pi/4 * cos(pi * x/(2*h))/h, esp))
  return(kh)
}

