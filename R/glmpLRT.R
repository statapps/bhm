#### penalized likelihood ratio test for glm

glmpLRT = function(x, ...) UseMethod("glmpLRT") 

glmpLRT.default=function(x, y, family, control, ...) {
  c0 = control$c0
  if (c0 == 0.5) {
    a = 1
    b = 1
    gc = 8
  } else if (c0 == .75) {
    a = 3
    b = 1
    gc = 64/3
  } else if (c0 == 0.25) {
    a = 1
    b = 3 
    gc = 64/3
  } else stop("c0 must be one of 0.25, 0.5, 0.75\n")
  control$ab = c(a, b, gc)
  
  z1 = NULL
  z2 = NULL
  
  p = ncol(x)
  #if (n.col > 4) cat("\nWarning: X has more than 4 columns, only the first two 
  #                   covariate will be used.")
  #if(n.col < 4) stop("Use glmpLRT(y~biomarker+trt+biomarker*trt)")
  if(p == 2) control$p1 = 0
  if(p > 2) z1 = x[, 3:p]
  p1 = control$p1

  # inverse cdf transformation of biomarker
  w = x[, 2]
  n = length(w)
  # Check for possible wrong order of formula
  if (length(unique(w)) < 4)
    stop("Either the formula is not in correct order of 'y~biomarker + trt' or the biomarker is not a continuous variable.")
  ew = ecdf(w)(w)

  if (control$p1 > 0) z2 = x[, c(1, 3:(3+p1-1))] ### with the interaction term
  else z2 = matrix(1, n, 1)   ### without the interaction term
  
  control$family = family
  #print(head(z1), digits = 3)
  #print(head(z2), digits = 3)
  # w: biomarker, z1: covariate, z2: interaction with w
  fit = .pLRTest(ew, y, z1, z2, control)
  fit$w = w
  fit$control = control
  
  if (control$method == "Bootstrap") {
    bLRT = .bglmTest(ew, y, z1, z2, control)
    fit$bpv = bLRT$bpv
    fit$bLR = bLRT$LR
  }

  fit$call = match.call()
  class(fit) = "glmpLRT"
  return(fit)
}

glmpLRT.formula = function(formula, family = binomial, data=list(...), 
      lambda = 15, c0 = 0.5, p1 = 1, method = c("pLRT", "Davies", "Bootstrap"), 
      B=10, K = 50, epsilon = 0.025,...){
  mf = model.frame(formula=formula, data=data)
  
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  method = match.arg(method)
  
  cq = seq(0.05, 0.95, epsilon)
  control = list(lambda = lambda, c0 = c0, p1 = p1, B = B, K = K, 
          cq = cq, method = method)

  fit = glmpLRT.default(x, y, family, control)
  fit$call = match.call()
  return(fit)
}

###### x: biomarker, y: outcome, z: treatment
.pLRTest = function(x, y, z1, z2, control) {
  n = length(x)
  
  a = control$ab[1]
  b = control$ab[2]
  gc = control$ab[3]
  c0 = control$c0              ## evaluate under the null H0: c = c0
  lambda = control$lambda
  family = control$family
  K = control$K
  
  ############ Davies test ######
  cq = control$cq #seq(0.1, 0.9, 0.05)
  J = length(cq)
  Mn = 0
  Rn = 0 
  
  Rx = rep(0, J)
  Mx = Rx
  c.max = NA  ### replace c.max = 0.5
  
  ### Null model without biomarker
  if(is.null(z1)) g0=glm(y~1,family=family)
  else g0=glm(y~z1,family=family)
  ml0 = logLik(g0)
  
  mpv = NULL
  sdf = NULL
  zeta = NULL
  #### Loop the profilelikelihood
  for (i in 1:J) {
    cx = cq[i]
    #w=exp(K*(x-cx))/(1+exp(K*(x-cx)))
    w = 1/(1+exp(-K*(x-cx)))   ### faster
    z2w = z2*w
    if(is.null(z1)) gm=glm(y~z2w, family=family)
    else gm=glm(y~z1+z2w, family=family)
    
    ml=logLik(gm)
    
    Ri = 2*(ml - ml0)
    Mi = Ri + 2*lambda*log((cx/c0)^a*((1-cx)/(1-c0))^b)
    Rx[i] = Ri
    Mx[i] = Mi
    if(Rn<Ri) {
      lp = predict(gm)
      gmx = gm
      c.max = cx
      Rn = Ri
      bn = gm$coef ##### coeficient to be used in bootstrap method
    }
    if(Mn<Mi) Mn = Mi
  }

  cLRT = cbind(cq, Mx, Rx)
  varNames = names(bn)
  if(control$method == 'pLRT') {
    #w0=exp(K*(x-c0))/(1+exp(K*(x-c0)))
    w0 = 1/(1+exp(-K*(x-c0)))   ###faster
    z2w0 = z2*w0
    #X = cbind(rep(1, n), z, w0, z*w0)
    #x2 = X[, 3:4]
    
    X = cbind(rep(1, n), z1, z2w0)

    nc = ncol(X)
    p2 = ncol(z2)-1
    idx1 = 1:(nc-p2-1)
    idx2 = (nc-p2):nc

    x2 = X[, idx2]  ### X matrix with biomarker term
    #print(head(X), digits=3)

    xb = predict(g0)
    eb = exp(xb)
    xp = eb/(1+eb)
    b2 = xp*(1-xp)
    vb = b2*(w0*(1-w0))^2*K*K
    
    An = t(X)%*%diag(b2)%*%X/n
    #cat("An = \n")
    #print(An, digits = 3)

    A22 = An[idx2, idx2]
    A12 = An[idx1, idx2]
    A11 = An[idx1, idx1]

    Aj = t(A12)%*%solve(A11)%*%A12
    #print(Aj, digits = 3)
    UU3 = t(x2)%*%diag(vb)%*%x2/n
    A22s = A22-UU3/(gc*lambda)
    H = A22s-Aj
    J2 = chol(A22-Aj)
    Ij = J2%*%solve(H)%*%t(J2)
    sdf = sum(diag(Ij))
    egn = eigen(Ij)$values
    if(min(egn) < 0) stop("Error: Negative eigen value here.\n")

    mu1 = sum(egn)
    mu2 = 2*sum(egn^2)
    zeta = mu2/(2*mu1)
    sdf = 2*mu1^2/mu2
    mpv = (1 - pchisq(Mn/zeta, sdf))
    #cat(' zeta = ', zeta, 'df=', sdf)
  }
  
  ##### Find p-values
  R2 = sqrt(Rx)
  s = 2
  V = sum(abs(diff(R2)))
  rpv = (1 - pchisq(Rn, s)) + V*Rn^{(s-1)/2}*exp(-0.5*Rn)*2^{-0.5*s}/gamma(0.5*s)
  if(rpv>1) rpv = 1
  return(list(mpv=mpv, rpv=rpv, cLRT = cLRT, df = sdf, varNames=varNames,
         zeta = zeta, lglk0 = ml0, loglik = Rn, Mn = Mn, gmx = gmx,
         c.max = c.max, coefficients = bn, linear.predictors = lp))
}
###### bootstrap
###### x: biomarker, y: outcome, z: treatment
.bglmTest = function(x, y, z1, z2, control) {
  if(is.null(z1)) p0 = 1
  else p0 = length(z1[1, ])+1
  m0 = .pLRTest(x, y, z1, z2, control)
  b1 = m0$coefficients[1:p0]
  Rn = m0$loglik
  n  = length(y)

  xb = cbind(rep(1, n), z1)%*%b1
  eb = exp(xb)
  xp = eb/(1+eb)

  B = control$B
  LR = rep(0, B)
  for (i in 1:B) {
    y0 = rbinom(n, 1, xp)
    fb = .pLRTest(x, y0, z1, z2, control)
    LR[i] = fb$loglik
  }
  bpv = mean(LR>Rn)
  return(list(LR=LR, bpv = bpv)) 
  #qtl = c(quantile(LR, 0.95), Rn)))
}

print.glmpLRT = function(x, ...) {
  cat("Call:\n")
  print(x$call)
  varNames = x$varNames
  p = length(varNames)
  p1 = x$control$p1
  cat("\nMain effect: ")
  cat(varNames[1:(p-p1-1)])
  cat("\nInteraction:", varNames[(p-p1):p])

  df  = x$df
  rpv = round(x$rpv*10000)/10000
  mpv = x$mpv
  bpv = x$bpv
  cat('\n\nDavis test p-value =', rpv, 'df =', p1+1, '\n')
  if(!is.null(mpv)) cat('pLRT  test p-value =', round(mpv*10000)/10000, 'df =', df, '\n')
  if(!is.null(bpv)) cat('Bootstrap test p-value =', x$bpv, '\n')
}

plot.glmpLRT = function(x, ..., scale = c('original', 'transformed')) {
  scale = match.arg(scale)
  wx = x$w
  lg0 = x$lglk0 
  cq = x$cLRT[, 1]
  Mn = x$cLRT[, 2] + lg0
  Rn = x$cLRT[, 3] + lg0

  w = quantile(wx, cq)
  w0 = switch(scale, original = w, transformed = cq)
  xlb = switch(scale, original = 'Biomarker', transformed = 'Biomarker Percentile')
  plot(w0, Rn, type = 'n', ylim = c(min(c(Mn, Rn)), max(Rn)), 
       xlab = xlb, 
       ylab = 'Log likelihood')
  lines(w0, Rn, lty = 1)
  lines(w0, Mn, lty = 2)
}
#plot(fit3)

summary.glmpLRT = function(object, ...) {
  if(object$rpv > 0.05) cat('There is not evidence to support the biomarker threshold and/or interaction effects')
  if(object$rpv <=0.05) {
    cat("The threshold parameter is estimated by\n")
    print(quantile(object$w, object$c.max))
    print(summary(object$gmx))
  }
}
