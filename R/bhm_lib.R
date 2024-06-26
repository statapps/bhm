### approximate a function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

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
    #### Please ensure the survival time data is sorted, run unext line to check for error
    #if(y[1, 1]<max(y[1, ])) stop("Survival time must be sorted from max to min.")
    # Do not run above line to speed up the algorithm.
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

gendat.glm =  function(n, c0, beta){
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

  kh = switch(kernel, gaussian = ifelse(ax < 5*h, dnorm(x, sd = h), esp),
    rectangular = ifelse(ax < h, 0.5/h, esp), 
    triangular = ifelse(ax < h, (1 - ax/h)/h, esp),
    epanechnikov = ifelse(ax < h, 3/4 * (1 - (ax/h)^2)/h, esp),
    biweight = ifelse(ax < h, 15/16 * (1 - (ax/h)^2)^2/h, esp),
    cosine = ifelse(ax < h, (1 + cos(pi * x/h))/(2*h), esp),
    optcosine = ifelse(ax < h, pi/4 * cos(pi * x/(2*h))/h, esp))
  return(kh)
}

### score function for the log likelihood using numerical derivative
numScore = function(func, theta, h = 0.0001, ...) {
  p = length(theta)
  score = rep(NA, p)
  for(i in 1:p) {
    theta1 = theta; theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    score[i] = (func(theta2, ...) - func(theta1, ...))/(2*h)
  }
  return(score)
}

### numerical Jacobian function 
numJacobian = function (func, theta, m, h = 0.0001, ...) {
  p = length(theta)
  jacobian = matrix(NA, m, p)
  for (i in 1:p) {
    theta1 = theta
    theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    jacobian[ ,i] = (func(theta2, ...) - func(theta1, ...))/(2*h)
  }
  return(jacobian)
}

### Hessian matrix for the log pseudo likelihood using numerical derivative
numHessian = function(func, theta, h = 0.0001, method = c("fast", "easy"), ...) {
  p = length(theta)
  H = matrix(NA, p, p)
  method = match.arg(method)

  ### easy to understand, evaluate func(...) [4*p*p] times
  if(method == "easy") {
    cat("Easy num hessian method\n")
    sc = function(theta, ...) {
      numScore(func = func, theta = theta, h = h, ...)
    }
    for(i in 1:p) {
      theta1 = theta; theta2 = theta
      theta1[i] = theta[i] - h
      theta2[i] = theta[i] + h
      H[i, ] = (sc(theta2, ...) - sc(theta1, ...))/(2*h)
    }
  }
  ### faster, evaluate func(theta, ...) [1 + p + p*p] times
  if (method == "fast") {
    cat("Fast num hessian method\n")
    f0 = func(theta, ...)
    for(i in 1:p) {
      theta1 = theta; theta2 = theta
      theta1[i] = theta[i] - h
      theta2[i] = theta[i] + h
      H[i, i] = (func(theta2, ...) - 2*f0 + func(theta1, ...))/h^2
    }
    for(i in 1:(p-1)){
      for(j in (i+1):p) {
        theta1 = theta; theta2 = theta
        theta1[i] = theta[i] - h; theta1[j] = theta[j] - h
        theta2[i] = theta[i] + h; theta2[j] = theta[j] + h
        H[i, j] = (func(theta2, ...) - 2*f0 + func(theta1, ...))/(2*h^2) - (H[i, i] + H[j, j])/2
        H[j, i] = H[i, j]
      }
    }
  }
  return(H)
}

### find roots for the multiple non-linear equations.
multiRoot = function(func, theta,..., verbose = FALSE, maxIter = 50, 
        thetaUp = NULL, thetaLow = NULL,
        tol = .Machine$double.eps^0.25) {
  alpha = 0.0001
  rho = 0.5
  U = func(theta, ...)
  m = length(U)
  p = length(theta)
  if (m == p) mp = TRUE  ### for m = p
  else mp = FALSE
  convergence = 0
  mU1 = sum(U^2)
  i = 1

  while(i < maxIter) {
    J = numJacobian(func, theta, m=m, ...)
    if(mp) dtheta = solve(J, U)
    else {
      tJ = t(J)       ## m*1-(p*m) x (m*p) x (p*m) x (m*1)
      dtheta = solve(tJ%*%J, tJ%*%U)
    }
    theta0 = theta
    lambda = 1
    delta = 1
    Ud = sum(U*dtheta)
    ########Linear search
    while (delta > 0) {
      theta = theta0 - lambda*dtheta
      if(!is.null(thetaUp)) theta = ifelse(theta > thetaUp, thetaUp, theta)
      if(!is.null(thetaLow)) theta = ifelse(theta < thetaLow, thetaLow, theta)
      U = func(theta, ...)
      mU = sum(U^2)
      delta = mU-mU1-alpha*lambda*Ud
      lambda = rho*lambda
      #if(verbose) cat('delta = ', delta, '\n')
    }
    dU = abs(mU1 - mU)
    mU1 = mU
    i = i + 1
    if(verbose) cat("||U|| = ", mU, "dU = ", dU, '\n')
    if ((mU < tol) | (dU < tol)) {
      convergence = 1
      if(mU > tol) convergence = 0
      break
    }
  }
  return(list(root = theta, f.root = U, iter = i, convergence = convergence))
}

#reverse rcumsum can be avoided if scored largest time to smallest time
#rcumsum=function(x) rev(cumsum(rev(x))) # sum from last to first

coxScoreHess = function(eb, delta, X, hess = FALSE) {
  ### eb = exp(x%*%beta)
  ### delta shall be sorted from smallest to largest.
  S0 = cumsum(eb)
  S1 = apply(eb*X, 2, cumsum)
  SX = delta * (X - S1/S0)
  score = colSums(SX)
  if(!hess) return(score)

  ### Sigma = Var(Score)
  Sigma = t(SX)%*%SX

  n = length(delta)
  p = ncol(X)
  SS1 = array(apply(S1, 1, function(x){return(x%*%t(x))}),c(p, p, n)) # ((p*p)*n)
  SS1 = aperm(array(SS1, c(p, p, n)), c(3, 1, 2))                     # (n*(p*p))

  Xt = apply(X, 1, function(x){return(x%*%t(x))})                    # X*t(X)
  X2 = array(Xt, c(p, p, n))
  ## multiply each X2(p, p, i) with eb[i],
  ## by change eb to a p*p*n array with each of ith pxp matrix = eb[i]
  X2eb = X2 * array(rep(eb, each = p*p), c(p, p, n))

  ## Sm is a upper triangular matrix of 1
  Sm = matrix(1, n, n)
  #Sm[lower.tri(Sm)] = 0
  Sm[upper.tri(Sm)] = 0

  ## calculate S2, a n*p*p array
  S2 = apply(X2eb, c(1, 2), function(x, y){return(y%*%x)}, Sm)
  H = colSums(delta*(S2/c(S0)-SS1/c(S0)^2), dims = 1)
  return(list(score = score, Sigma = Sigma, H = H))
}
