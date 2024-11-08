############## adaptive reject sampling ########################################
### ars is an adaptive reject sample method that make use of the 
### log likelihood function directly, It can choose the envelope 
### function automatically and it is not necessary to find 
### the max.h as in the other reject sample function
###
### this method allows use the sample directly from the log-likelihood function 
### with an unknown normalized constant. 
###
### Let h(x) = log(f(x)) defines the log-likelihood function
###
##################### main ars function ########################################
ars = function(h, ..., lower = -5, upper = 5, sigma.offset = 0.05, 
               mu.offset = 0.03, verbose = FALSE) {
  #print(lower)
  #print(upper)
  width = (upper - lower)/50
  #print(width)
  x0  = seq(lower, upper, width)
  
  hx0 = h(x0, ...)
  #print(h)
  max.h = max(hx0)
  hx0 = hx0 - max.h
  #print(hx0)
  mu0 = x0[hx0==0]
  
  lbd = -15
  xg  =  x0[hx0>lbd]
  hg  = hx0[hx0>lbd]
  lower = min(xg)
  upper = max(xg)

  eps = 1e-6
  ### numerical second order derivative of the log likelihood = 1/sigma^2
  ### inflated sigma by sigma.offset
  sigma = sqrt(abs(eps*eps/(2*max.h - h(mu0+eps, ...)-h(mu0-eps, ...)))) + sigma.offset
  pK = - dnorm(0, 0, sigma, log = TRUE) + mu.offset

  accept = FALSE
  ### generate N normal samples
  while(!accept) {
    x  = rnorm(1, mean = mu0, sd = sigma) #candidate distribution
    if(x < lower) x = lower
    if(x > upper) x = upper
    dx = dnorm(x, mean = mu0, sd = sigma, log=TRUE) + pK #envelope
    hx = h(x, ...) - max.h

    alpha = exp(hx-dx)
    if (runif(1) < alpha) accept = TRUE
  }
  if(verbose) {
    ng = dnorm(xg, mean = mu0, sd = sigma, log = TRUE) + pK
    plot( xg, exp(ng), type = 'l', ylim=c(0, 1), xlab = 'x', ylab = 'density')
    lines(xg, exp(hg), lty = 2)
    
    cat("Possible range:", lower, "to", upper, 
        "\nProposal mu =", mu0, "sigma =", sigma)  
    cat(". x sample = ", x, "\n")
  }
  return(x)
}
