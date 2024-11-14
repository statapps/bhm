############## adaptive reject sampling ########################################
### ars and arns adaptive reject sample method that make use of the 
### log likelihood function directly, It can choose the envelope 
### function automatically.
###
### this method allows use the sample directly from the log-likelihood function 
### with an unknown normalized constant. 
###
### Let h(x) = log(f(x)) defines the log-likelihood function
###
##################### main ars function ########################################
ars = function(logpdf, n = 1, lower = -14, upper = 15, x0 = 0, ...) {
  ### the log pdf can not be too small
  xmin = log(.Machine$double.xmin)
  eps = 1e-7
  samples = rep(NA, n)

  cx = c(lower, x0, upper)
  h  = logpdf(cx, ...)
  max.h = max(h)

  ### redefine logpdf function 
  fn = function(x) logpdf(x, ...) - max.h

  while(fn(lower) < xmin) lower = lower + 0.1
  while(fn(upper) < xmin) upper = upper - 0.1

  Tk = unique(c(lower, x0, upper))
  K  = length(Tk)
  h  = fn(Tk)
  dh = (fn(Tk+eps)-fn(Tk))/eps

  if(dh[1] < 0) cat('lower = ', lower, 'dh(lower) shall be positive\n')
  if(dh[K] > 0) cat('upper = ', upper, 'dh(upper) shall be negative\n')

  i = 1
  while(is.na(samples[n])) {
    ### proposed a sample xp
    #print(rbind(Tk, h, dh))
    xp = rpwexp(h, dh, Tk, lower, upper)

    ### evaluate the value of envlop function at xp
    ### at function gu(x) = h + (x - Tk)*dh
    ux = min(h + (xp - Tk)*dh) 
    gx = fn(xp)
    #cat('xp = ', xp, 'gx = ', gx, 'ux = ', ux, '\n')
    if(runif(1) < exp(gx-ux)) {
      samples[i] = xp
      i = i + 1
    } else {
      Tk = sort(c(xp, Tk))
      h  = fn(Tk)
      dh = (fn(Tk+eps)-fn(Tk))/eps
    }
  }
  return(samples)
}

### sample a piecewise exponential dist with lines go through (Tk, h(Tk))
### slope dh, intersect at z and intercept at c0:
### u(x) = h + (x - Tk)*dh, 
### z  = -(diff(h)-diff(Tk*dh))/diff(dh)
### c0 = h - Tk*dh
rpwexp = function(h, dh, Tk, lower, upper) {
  z  = c(lower, -(diff(h)-diff(Tk*dh))/diff(dh), upper)
  c0 = h - Tk*dh
  K = length(h) 
  ## int_a^b exp(c0 + dh*x)dx = (exp(c0+dh*b) - exp(c0+dh*a))/dh
  f = (exp(c0 + dh*z[2:(K+1)]) - exp(c0 + dh*z[1:K]))/dh  
  f = f/sum(f)
  j = sample(1:K, 1, prob = f)
  dj= dh[j]
  u = runif(1)
  x = log(u*exp(dj*z[j+1]) + (1-u)*exp(dj*z[j]))/dj
  return(x)
}

### adaptive reject normal sampling, using a normal distribution as an envelope function
arns = function (logpdf, n = 1, lower = -5, upper = 5, sigma.offset = 0.05,
          fx.offset = 0.03, K = 100, verbose = FALSE, ...) {

  width = (upper - lower)/K
  x0 = seq(lower, upper, width)
  hx0 = logpdf(x0, ...)
  max.h = max(hx0)
  hx0 = hx0 - max.h

  h = function(x) logpdf(x, ...) - max.h
  mu0 = x0[hx0==0]
  lbd = -15
  xg  =  x0[hx0>lbd]
  hg  = hx0[hx0>lbd]
  lower = min(xg)
  upper = max(xg)

  eps = 1e-7
  ### numerical second order derivative of the log likelihood = 1/sigma^2
  ### inflated sigma by sigma.offset
  # sigma=sqrt(abs(eps*eps/(2*max.h - h(mu0+eps, ...)-h(mu0-eps, ...)))) + sigma.offset
  sigma = sqrt(abs(eps*eps/(-h(mu0+eps) - h(mu0-eps)))) + sigma.offset

  pK = - dnorm(0, 0, sigma, log = TRUE) + fx.offset

  samples = rep(NA, n)
  i = 1
  ### generate N normal samples
  while (is.na(samples[n])) {
    x  = rnorm(1, mean = mu0, sd = sigma) #candidate distribution
    if(x < lower) x = lower
    if(x > upper) x = upper
    dx = dnorm(x, mean = mu0, sd = sigma, log=TRUE) + pK #envelope
    hx = h(x)

    alpha = exp(hx-dx)
    if (runif(1) < alpha) {
      accept = TRUE
      samples[i] = x
      i = i + 1
    }
  }
  if(verbose) {
    ng = dnorm(xg, mean = mu0, sd = sigma, log = TRUE) + pK
    plot( xg, exp(ng), type = 'l', ylim=c(0, 1), xlab = 'x', ylab = 'density')
    lines(xg, exp(hg), lty = 2)
    
    cat("Possible range:", lower, "to", upper, 
        "\nProposal mu =", mu0, "sigma =", sigma)  
    cat(". x sample = ", x, "\n")
  }
  return(samples)
}
### end of functions 
