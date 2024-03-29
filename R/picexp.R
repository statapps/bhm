##########Note: git pull before working on the code.
hpicexp = function(x, rate, cuts, index=NULL) {
  if(is.null(index)) index = findInterval(x, cuts)
  return(rate[index])
}

### derivative of h(t) w.r.t baseline rate parameter.
.dhpicexp = function(x, rate, cuts) {
  ### p is the col number of rate, which is n.cuts - 1
  p = length(cuts)-1
  n = length(x)
  
  ### dh_i = 1 if x is in the ith interval, 0 otherwise
  dh = matrix(0, p, n)
  
  ### access element of dh matrix as a vector
  s = seq.int(0, n-1, 1)*p
  index = findInterval(x, cuts) + s
  dh[index] = 1
  return(t(dh))
}

Hpicexp = function(x, rate, cuts, index=NULL) {
  if(is.null(index)) index = findInterval(x, cuts)
  dx = diff(cuts)*rate
  cdx = c(0, cumsum(dx))
  ### cumlative hazard up to last cup point
  Ha = cdx[index]
  
  ### rest of cum haz 
  Hb = (x-cuts[index]) * rate[index]
  H = Ha+Hb
  return(H)
}

.dHpicexp = function(x, rate, cuts) {
  ### p is the col number of rate, which is n.cuts - 1
  p = length(cuts)-1
  n = length(x)
  
  dH = matrix(0, n, p)
  for(j in 1:p) {
    a = ifelse(x < cuts[j+1], x, cuts[j+1]) - cuts[j]
    dH[, j] = ifelse(a>0, a, 0)
  }
  return(dH)
}

ppicexp = function(q, rate=1, cuts=c(0, 10), lower.tail = TRUE, index = NULL) {
  H = Hpicexp(q, rate, cuts, index)
  s = exp(-H)
  ## If lower.tail then return F = P(T <= q)
  if(lower.tail) return (1-s)
  else return (s)
}

qpicexp = function(p, rate=1, cuts=c(0, 10), lower.tail=TRUE) {
  if(lower.tail==FALSE) p = (1-p)
  qc = cuts
  cuts.n = length(cuts)
  qc[cuts.n] = cuts[cuts.n] - 1E-9
  Hc = Hpicexp(qc, rate, cuts)
  lp = -log(1-p)
  index = as.numeric(cut(lp, Hc))
  qa = cuts[index]
  
  qb = (lp-Hc[index])/rate[index]
  return(qa+qb)
}

rpicexp = function(n, rate = 1, cuts=c(0, 10)) {
  return(qpicexp(runif(n, 0, 1), rate, cuts))
}

dpicexp = function(x, rate=1, cuts=c(0, 10), log = FALSE) {
  index = findInterval(x, cuts)
  h = hpicexp(x, rate, cuts, index)
  s = ppicexp(x, rate, cuts, lower.tail=FALSE, index)
  if(log) return(log(h) + log(s))
  else return(h*s);
}

############ Piecewise exponential regression ###########
### log likelihood function 
.logPic = function(lambda, y, cuts) {
  time = y[, 1]
  event = y[, 2]
  
  h = hpicexp(time, lambda, cuts)
  H = Hpicexp(time, lambda, cuts)
  ell = event*log(h) - H
  return(-sum(ell))
}

.scorePic = function(lambda, y, cuts) {
  time = y[, 1]
  event = y[, 2]
  
  h = hpicexp(time, lambda, cuts)
  dh = .dhpicexp(time, lambda, cuts)
  dH = .dHpicexp(time, lambda, cuts)
  dl = -(event/h*dh - dH)
  return((apply(dl, 2, sum)))
}

picreg = function(...) {
  UseMethod("picreg")
}

picfit = function(y, cuts=c(0, 10)) {
  if(cuts[1] !=0) stop("cuts[1] must be 0.")
  if(min(0, diff(cuts))<0) stop("custs cannot decrease.")
  time = y[, 1]
  event = y[, 2]
  p0 = length(cuts)-1
  lambda = rep(0.2, p0)
  
  ml = optim(lambda, .logPic, .scorePic,lower = 1e-9, method = "L", hessian = TRUE,
             y = y, cuts=cuts)
  lambda = ml$par
  var = solve(ml$hessian)
  sd = sqrt(diag(var))
  fit = list(coefficients=ml$par, var = var, sd = sd, 
                    logLik=-ml$value)
  class(fit) = c("picreg")
  return(fit)
}

print.picreg=function(x, digits=3,...) {
  cat("picreg: AFT model with piecewise exponential\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat('Coefficients and baseline parameter\n')
  bc = cbind(unname(t(x$coefficients)))
  colnames(bc) = t(x$varNames)
  print(bc, digits=digits)
  cat("log likelihood = ", x$logLik, "\n")
}

summary.picreg = function(object, alpha = 0.05,...){
  zscore = qnorm(1-alpha/2)
  sd = object$sd
  theta = t(unname(t(object$coefficients)))
  
  qtl = cbind(theta-zscore*sd, theta+zscore*sd)
  
  TAB1 = cbind(theta, sd, theta/sd, qtl[,1], qtl[,2], 2*pnorm(-abs(theta/sd)))
  colnames(TAB1) = c("Estimate", "Std.Err", "Z value", "95% CI(Low)", "95% CI(Up)", "Pr(>z)")
  rownames(TAB1) = object$varNames
  
  #lp = -(object$linear.predictors)
  #cidx = concordance(object$y~lp)
  #iter = object$iter
  results = list(call=object$call,TAB1=TAB1, logLik=object$logLik) #cidx = cidx)
  class(results) = "summary.picreg"
  return(results)
}

### derivative of h(t) w.r.t time parameter.
### approximated by linear
.dhtpicexp = function(x, rate, cuts) {
  h = 0.05
  cut1 = cuts+h
  cut2 = sort(c(cuts, cut1))
  index = findInterval(x, cut2)
  drate = diff(c(0, rate))/h
  dht = ifelse(index%%2, drate[ceiling(index/2)], 0)
}

#########Useless
#.qpicexp.slow = function(p, rate=1, cuts=c(0, 10)) {
#  p = as.matrix(p)
#  max.cut = max(cuts) - 1E-9
#  if(max(p) > ppicexp(max.cut, rate, cuts)) 
#    stop("Error: p too large, consider a large value for the max cut point.")
#  f = function(x, p0) {(ppicexp(x, rate, cuts) - p0)}
#  r = function(p0) uniroot(f, c(0, max.cut), p0 = p0)$root
#  q = apply(p, 1, r)
#  return(q)
#}

#picexp = list(name = "Picexp", 
#        init = function(x, weights, ...) {
#          c(mean(x), var(x))},
#        density = function(x, parms) {
#       }
#        deviance= function(...) stop('deviance residuals not defined')
#      )
