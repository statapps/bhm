#generate regression coefficients
thm.fit <-function(x, y, family, cx){
  #Remove the intercept term
  c.n = length(cx)
  n.col = length(x[1, ])
  x = x[, -1]
  if(n.col == 2) x = matrix(x, ncol = 1)
  w = x[, 1]
  if (c.n == 1) x[, 1] = ifelse(cx <= w, 1, 0)
  if (c.n == 2) x[, 1] = ifelse((cx[1]<=w) & (w<=cx[2]), 1, 0)
  
  if (length(grep(":", colnames(x)[3])) == 1) 
    x[, 3]=x[, 1]*x[, 2]
  else if(length(grep(":", colnames(x)[2])) == 1) 
    x[, 2] = x[, 1]*x[, 2]

  x.c_ = x
  if (family == "surv") {
    fit = coxph(y~x.c_)
    fit$converged = TRUE
  }
  else 
    fit = glm(y~x.c_, family=family)
  
  return(fit)
}

#calculate the joint probability (log likelihood) of y
thm.lik = function(x, y, family, beta, q, cx, control){
  c.n = control$c.n
  beta0 = control$beta0
  sigma0.inv = (control$sigma0.inv)
  n.col = length(x[1, ])

  w = x[, 2]
  if (c.n == 1) x[, 2] = ifelse(cx <= w, 1, 0)
  if (c.n == 2) x[, 2] = ifelse((cx[1]<=w) & (w<=cx[2]), 1, 0)

  if(length(grep(":", colnames(x)[4])) == 1) 
    x[, 4] = x[, 2]*x[, 3]
  else if(length(grep(":", colnames(x)[3])) == 1)
    x[, 2] = x[, 2]*x[, 3]

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
  lik = lik + lik2 - 0.5*t(beta-beta0)%*%sigma0.inv%*%(beta-beta0)
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
  h0 = 1
  h = h0*exp(X%*%beta)
  stime = rexp(n, h)               #Failure time.
  endstudy = runif(n, 2, 4)
  cstatus = ifelse(stime>endstudy, 0, 1) ##Censored at end of study time.
  stime = ifelse(stime>endstudy, endstudy, stime)
  dat = cbind(stime, cstatus, z, x)
  return(dat)
}

glm.gendat =  function(n, c0, beta){
  n1 = n/2
  x = runif(n, 0, 1)
  c.n = length(c0)
  
  if (c.n == 1) x1 = ifelse(x > c0, 1, 0)
  if (c.n == 2) x1 = ifelse((c0[1] < x) & (x < c0[2]), 1, 0)

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
    p[i] = sum(x<x[i])
  }
  p = (p-0.5)/n
  return(p)
}
