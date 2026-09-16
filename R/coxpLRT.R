#source("library.R")
#library(matrixStats)

coxpLRT = function(XZ, ...) UseMethod("coxpLRT")
### mmList = model matrix list of covariate interacts with x
### x is the bio-marker (x1, x2, x3)
### z is the the rest of covariates in the model, 
###   such as the covariate main effect, age, state , etc. 
###   formulaList is a vector of covariate interacts with x1, x2,  x3 ...
coxpLRT.formula = function(formula, formula2 = NULL, formulaList = NULL, 
          data=list(...), c0 = 0.5, lambda = 25, K = 50, x = FALSE, 
          bootstrap = 0, ...) {
  if(x) returnX = TRUE else returnX = FALSE

  mf = model.frame(formula=formula, data=data)
  y = model.response(mf)  # survival outcomes including censoring indicator.

  ## check survival object
  if (is(y, "Surv")) {
    idx = order(y[, 1], decreasing = TRUE)
    y = y[idx, ]
    data = data[idx, ] ### sort the data
  } else stop("y shall be a Surv object.")
  
  XZ  = .getXZ(formula, formula2, formulaList, data, c0, K)
  fit = coxpLRT.default(XZ, y, lambda = lambda)
  if(bootstrap>0) {
    fit$boot = .resbootCox(XZ, y, bootstrap)
  } else {fit$boot$pValue = 1}
  return(fit)
}

####### penalized Likelihood ratio test for the Cox PH model 
####### a good initial value of theta help speed up the algorithm
coxpLRT.default = function(XZ, y, lambda = 25, naive.test = FALSE, ...) {
### theta:1:p0 for alpha * x, 
###       p0+      1 to p0+pm[1] for beta1, 
###       p0+pm[1]+1 to pm+pm[2] for beta2,
###       p0+pm[2]+1 to pm+pm[3] for beta3 etc.  
### y     survival object
### 
### We need x, z, z2 as matrix
  K  = XZ$K
  x  = as.matrix(XZ$x)    ### x    : n x m biomarker matrix
  z  = as.matrix(XZ$z)    ### z    : covariate main effect
  z2 = as.matrix(XZ$z2)   ### z2   : covariate interact with biomarkers, intercept as biomarker main effect
  pm = XZ$pm              ### pm   : pm[i-1]+1 to pm[i] is the columns index of z2 that
                          ###        interacts with biomarker
  #cxterm = .plterm(cx0)  ### Calculate the penalty term, to be used in Gc below
  cxterm = XZ$cxterm
  
  n  = nrow(x)
  m  = ncol(x)
  p2 = pm[m+1] ######degree of freedom for the Wald and robust test
  if(length(pm) !=(m+1)) stop("pm must be a vector of length m+1")
  p0 = ncol(z)  ### z must be a matrix
  p  = p0 + p2

  ### Null model without the bio-marker
  if (is.null(z)) {
    g0=coxph(y~1)
  } else {
    g0=coxph(y~z)    # Covariates (other than biomarker) main effects. Model under H0.
  }
  
  # theta[1:p0]=alpha, theta[(p0+1) : p] = beta, X'*alpha + Zc'*beta
  # theta[(p+1):(p+m)] = cut points
  cx0   = XZ$c0
  theta = c(rep(0, p), cx0)
  theta[1:p0] = g0$coef
  
  tm = y[, 1]                #########Survival time shall be sorted
  if((tm[1]<tm[2])|(tm[2]<tm[3])) stop("Survival time is not sorted: T1>T2>T3...>Tn.")
  delta = y[, 2]
  
  thetaUp  = c(rep( 10, p), rep(0.95, m))
  thetaLow = c(rep(-10, p), rep(0.05, m))
  
  obj = optim(theta, .coxpLoglik, gr = .coxpLRTScore, 
                event = delta, x=x, z=z, z2=z2, lambda = lambda, 
                pm = pm, K=K, cx = cxterm, 
                method = "L-BFGS-B", 
                lower = thetaLow, upper = thetaUp)
  theta = obj$par

  ### plot profile likelihood
  #cplikPlot(theta, x, delta, z, z2, lambda, pm, K, cxterm)

  ### cxterm will be used in the calculation of pLRT, 
  ### cutpoint is obtained from the MLE of theta 
  sc = .coxLoglikScore(theta, delta, x, z, z2, lambda, pm, K, cxterm, score=TRUE)
  
  ml0 = logLik(g0)          ### ml0 is the logLik under H0
  eb = exp(predict(g0))
  
  ml = sc$logLik            ### ml  is the MLE under H1 
  Z_ = sc$Zc                ### new Z_, with optimal cutpoint
  varNames  = XZ$varNames
  colnames(Z_) = varNames
  gmx   = coxph(y~Z_)
  c.max = theta[(p+1):(p+m)] ### lambda = 0 return c.max

  ### return model without penalty term, 
  ### this does not control the overall type I error
  if(lambda == 0) {
    lr1 = logLik(gmx)
    lr  = 2*as.numeric(lr1-ml0)

    ### LRT Wald and Robust tests
    if (naive.test) {
      df1 = attr(lr1, 'df')
      df0 = attr(ml0, 'df')
      dfr = df1-df0
      sel = (p0+1):p
      V = vcov(gmx)
      b = coef(gmx)[sel]
      xk2 = t(b)%*%solve(V[sel, sel])%*%b
      
      ### Robust rest
      gmxr = coxph(y~Z_, robust = TRUE)
      xkr = t(b)%*%solve(gmxr$var[sel, sel])%*%b
      #print(gmxr$naive.var)
      
      pv1 = pchisq(lr,  dfr, lower.tail = FALSE)
      pv2 = pchisq(xk2, dfr, lower.tail = FALSE)
      pv3 = pchisq(xkr, dfr, lower.tail = FALSE)
      #lrt = structure(lr, class='logLik', df = dfr, nobs = n)
      pValue = c(pv1, pv2, pv3)
    } else {
      pValue = c(NA, NA, NA)   # skip Wald/robust entirely for bootstrap speed
    }
    return(list(fit1 = gmx, fit0 = g0, lr = lr, pValue = pValue,
           tests = c("nLRT", "nWald", "nRobust"), 
           c.max = c.max, theta = c(gmx$coef, c.max), 
           z = z, Zc = Z_, 
      message = "Warning: With lambda = 0, this test does not control the Type I error.")
    )
  }
  
  bn = gmx$coef
  names(bn) = varNames
  lp = predict(gmx)
  ebx = exp(lp)
  
############## below is to calculate pLRT, Wald and robust test under H0 ###########
#### theta  is the MPPL estimate in the manuscript
#### theta1 is the \tilde{theta} in the manuscript
  Zc     = z
  z2dw0  = NULL
  ##### crate the model matrix Zc under H0 with cut point to c0
  for(i in 1:m) {
    zi    = z2[, (pm[i]+1):pm[i+1]]
    w0    = 1/(1+exp(-K*(x[, i]-cx0[i])))  ### faster
    dw0   = -K*w0*(1-w0)                   ### faster
    Zc    = cbind(Zc,    zi*w0)            #z2w0  = zi*w0
    z2dw0 = cbind(z2dw0, zi*dw0)           #z2dw0 = zi*dw0, we need this for U3 below
  }
  
  gc     = cxterm[4, ] ###For cx0 = 0.5 gc = 8. For cx0 = 0.25, 0.75, gc = 64/3 
  theta1 = theta
  theta1[1:p0]        = g0$coef
  theta1[(p0+1):p]    = 0
  theta1[(p+1):(p+m)] = cx0
  #print(theta1)

  #alpha = theta[1:p0]  #alpha is not need here
  beta  = theta[(p0+1):p]
  gma   = theta[1:p]  ### regression coefs

  Mn = 2*(ml - ml0)
  #print(Mn)
  
  mpv = NULL
  sdf = NULL
  zeta = NULL
  
  S0 = cumsum(eb)
  S1U3 = apply(z2dw0*eb, 2, cumsum) #n*p2
  U3  = colSums(delta*(z2dw0 - S1U3/S0))/sqrt(n)

  ## U3m is a p2 x m matrix
  U3m = matrix(0, p2, m)
  for(i in 1:m) {
    idx = (pm[i]+1):pm[i+1]
    U3m[idx, i] = U3[idx]
  }

  ### UU3 is U3*U3, a p2*p2 matrix from formula U3
  UU3 = U3m%*%t(U3m)     ## p2*p2.
  
  ### coxScoreHess is a function in the lpl package to calculate 
  ### score and Hess matrix for cox model
  cfit = coxScoreHess(X=Zc, y = y, exb = eb, hess=TRUE)
  An = cfit$H/n    ########  A_n = I_n is the Hessian under H0

  idx1 = 1:p0
  idx2 = (p0+1):p
  
  ### An = I_n(tilde theta) in the manuscript
  A22 = matrix(An[idx2, idx2], nrow=length(idx2))
  A12 = matrix(An[idx1, idx2], nrow=length(idx1))
  A11 = matrix(An[idx1, idx1], nrow=length(idx1))
  
  A11invA12 = solve(A11)%*%A12
  
  Aj = t(A12)%*%A11invA12
  #print(Aj, digits = 3)
  
  Jn = A22 - Aj
  Gc = diag(rep(1/(gc*lambda), diff(pm)))
  Qn = Jn - UU3%*%Gc   ## Jn - UU3/(8*lambda) for c0 = 0.5
  #print(UU3%*%Gc)
  
  #H = A22s-Aj
  #J2 = chol(A22-Aj)
  #Ij = J2%*%solve(H)%*%t(J2)
  
  J2 = chol(Jn)  ######## J2 = (Jn)^0.5, 
  Ij = J2%*%solve(Qn)%*%t(J2)
  
  sdf = sum(diag(Ij))
  egn = eigen(Ij)$values
  if(min(egn) < 0) stop("Error: Negative eigen value here, use a different lambda and K.\n")
  mu1 = sum(egn)
  mu2 = 2*sum(egn^2)
  zeta = mu2/(2*mu1)
  sdf = 2*mu1^2/mu2

  mpv = (1 - pchisq(Mn/zeta, sdf)) ## p-value for PLRT.
  #cat(' zeta = ', zeta, 'df=', sdf, "Eigen value =  ", egn, "pvalue = ", mpv, "\n")

  ######## Wald test
  Wn1 = n*t(beta)%*%Qn%*%solve(Jn)%*%Qn%*%beta
  
  ###recalculate under theta hat, ebx = exp(beta hat * Z_) for robust rest
  cfitx = coxScoreHess(X = Z_, y = y, exb = ebx, hess=TRUE)

  Sn = cfitx$Sigma/n
  S22 = matrix(Sn[idx2, idx2], nrow=length(idx2))
  S12 = matrix(Sn[idx1, idx2], nrow=length(idx1))
  S11 = matrix(Sn[idx1, idx1], nrow=length(idx1))
  
  Js = t(A11invA12)%*%S11%*%A11invA12 + S22 - 2*t(A11invA12)%*%S12
  Wn2 = n*t(beta)%*%Qn%*%solve(Js)%*%Qn%*%beta
  Wn  = c(Wn1, Wn2)
  wpv = 1 - pchisq(Wn, p2)
  test = c(Mn/zeta, Wn1, Wn2)
  
  fit= list(mpv=mpv, df = c(sdf, p2, p2), test = test, 
            xNames=XZ$xNames, pValue = c(mpv, wpv), 
            testName = c('pLRT', 'Wald', 'Robust'),
            zNames = XZ$zNames, varNames=varNames, pm = pm, 
            zeta = zeta, Zx = Z_, theta = theta,
            iter = obj$iter, convergence = obj$convergence, 
            lglk0 = ml0, loglik = ml, Mn = Mn, coxCmx= gmx,
            c.max = c.max, coefficients = bn, linear.predictors = lp
            )
  class(fit) = "coxpLRT"
  return(fit)
}

### Approximate function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }
.resbootCox = function(XZ, y, bootstrap) {
  lambda = 0
  n = nrow(y)
  
  ### turn of naive test to speed up
  LR = coxpLRT.default(XZ, y, lambda, naive.test = FALSE)
  mLRT = LR$lr
  z  = LR$z
  Zc = LR$Zc
  B = bootstrap
  
  cfit0 = LR$fit0
  cfit1 = LR$fit1
  
  S0 = basehaz(cfit1, centered = FALSE)
  
  beta  = cfit1$coeff
  beta0 = cfit0$coeff
  
  exb  = exp(Zc%*%beta)
  exb0 = exp(z%*%beta0)
  
  #estimated residuals under the alternative hypothesis
  chaz = S0$hazard
  tm0  = S0$time
  xchz = .appxf(chaz, tm0, y[, 1])
  uhat = (exp(-xchz))^exb
  
  i=1
  LRb = rep(0, B)
  XZb = XZ
  while (i<=B) {
    idx = sample(1:n, replace=TRUE)
    
    uStar = uhat[idx]
    yStar = y[idx, ]
    tmStar=1-uStar^(1/exb0)
    yStar[, 1] = tmStar
    
    ### sort the time Star of the bootstrap sample
    idx2   = order(tmStar, decreasing = TRUE)
    yStar  = yStar[idx2, ]
    XZb$x  = XZ$x[idx2,  ]
    XZb$z  = XZ$z[idx2,  ]
    XZb$z2 = XZ$z2[idx2, ]
    
    lStar = try(coxpLRT.default(XZb, yStar, lambda, naive.test = FALSE))
    if(is(lStar, "try-error")) next
    
    LRb[i] = lStar$lr
    i=i+1
  }
  pValue=sum(LRb>mLRT)/B
  
  return(list(pValue = pValue, mLRT = mLRT, LRb = LRb))
}

resbootCox = function(formula, formula2 = NULL, formulaList = NULL, 
            data=list(...), bootstrap = 10, ...) {
  fit = coxpLRT.formula(formula, formula2 = formula2, 
          formulaList = formulaList, data=data,bootstrap = bootstrap) 
  return(fit$boot)
}

#cplikPlot = function(theta, x, delta, z, z2, lambda, pm, K, cxterm) {
#  m  = ncol(x)
#  p0 = ncol(z) 
#  p  = p0 + pm[m+1]
#  cq = seq(-0.1, 0.9, 0.02);  lk = cq
  
#  for(i in 1:length(cq)) {
#    thetai = theta; thetai[2] = cq[i]
#    lk[i] = .coxLoglikScore(thetai, delta, x, z, z2, lambda, pm, K, cxterm)
#  }
#  plot(cq, 2*lk, type = 'l', xlab = 'beta1', ylab = '2logLik')
  #abline(v = 0, lty = 2)
#  b = theta[2]
#  c = 2*max(lk)
#  v = vcov(fit$coxCmx)[2, 2]
#  lk2 = c - n*v*(cq-b)^2
#  print(lk2)
#  lines(cq, lk2, lty = 2)
#}

plot.coxpLRT = function(x, ...) {}

print.coxpLRT = function (x, ...) {
  cat("Call:\n")
  print(x$call)

  #p = length(varNames)
  p1 = x$control$p1
  cat("\nMain effect: ")
  cat(x$zNames)
  
  cat("\nInteraction:", x$xNames, "\n")
  df = x$df
  rpv = round(x$rpv * 10000)/10000
  pv = x$pValue
  bpv = x$boot$pValue
  #Wn  = x$test[2]
  #dfw = x$df[2]
  #wpv = x$wpv
  if (!is.null(pv[1])) 
    cat("pLRT   test p-value =", round(pv[1] * 10000)/10000, 
        "df =", df[1], "\n")
  cat("Wald   test p-value =", pv[2], "df =", df[2], "\n")
  cat("Robust test p-value =", pv[3], "df =", df[2], "\n\n")
  if (!is.null(bpv)) 
    cat("Bootstrap test p-value =", x$bpv, "\n")
  cat("Cox PH model with the optimal cut point = ", x$c.max, " :\n")
  print(x$coxCmx)
}


###Process the value of penalty term using c0 for the pLRT
.plterm = function(c0) {
  m = length(c0)
  cxterm = matrix(1, 4, m)
  for(i in 1:m) {
    if (c0[i] == 0.5) {
      cxterm[1, i] = 1
      cxterm[2, i] = 1
      cxterm[3, i] = 0.5
      cxterm[4, i] = 8
    } else if (c0[i] == .75) {
      cxterm[1, i] = 3
      cxterm[2, i] = 1
      cxterm[3, i] = 0.75
      cxterm[4, i] = 64/3
    } else if (c0[i] == 0.25) {
      cxterm[1, i] = 1
      cxterm[2, i] = 3
      cxterm[3, i] = 0.25
      cxterm[4, i] = 64/3
    } else stop("c0 must be one of 0.25, 0.5, 0.75\n")
  }
  return(cxterm)
}

### obtain the biomarker x, convariate z and XZ interaction term from the data
### using the cutpoints cx
.getXZ = function(formula, formula2, formulaList, data, cx, K) {
  x  = .matX(formula, data)
  m  = ncol(x)
  pm = rep(0, m+1)
  z  = .matX(formula2, data)
  
  z2 = NULL  ### z2 is for the covariates interact with x
  for(i in 1:length(formulaList)) {
    z2i = model.matrix(formulaList[[i]], data = data)
    z2  = cbind(z2, z2i)
    pm[i+1] = ncol(z2i)
  }
  pm = cumsum(pm)
  p  = pm[m+1] + ncol(z)
  
  if(length(cx) == 1) cx = rep(cx, m)
  cxterm = .plterm(cx) ###Calculate the penalty term
  
  Zc = z
  
  ### var names for Zc 
  zNames = colnames(z)
  xNames = colnames(x)
  varNames= zNames
  
  for(i in 1:m) {
    zi  = z2[, (pm[i]+1):(pm[i+1])]     ### use cx as the cut point
    w0  = 1/(1+exp(-K*(x[, i]-cx[i])))  ### faster
    Zc  = cbind(Zc, zi*w0)
    
    dw0 = -K*w0*(1-w0)                  ### faster
    
    ziNames = colnames(zi)[-1]  ### remove the intercept name
    varNames = c(varNames, xNames[i])
    for(j in 1:length(ziNames)) varNames = c(varNames, paste(xNames[i], '*', ziNames[j], sep = ''))
  }
  return(list(x=x, z = z, z2 = z2, Zc = Zc, m = m, p = p, pm = pm, K = K,
              c0 = cx, cxterm = cxterm, 
              varNames = varNames, xNames=xNames, zNames=zNames))
}

### get the covariate X (without the Intercept) used for the coxph model
.matX = function(formula, data){
  mf = model.frame(formula=formula, data=data)
  x  = model.matrix(attr(mf, "terms"), data = mf)
  x  = x[, -1, drop = FALSE]
  return(x)
}

## return score of the Cox model, to be used in the optimization function
## use cache to speed up, see below
#.coxpLRTScore = function(theta, event, x, z, z2, lambda, pm, K, cx)
#  -.coxLoglikScore(theta, event, x, z, z2, lambda, pm, K, cx, score = TRUE)$score

#.coxpLoglik = function(theta, event, x, z, z2, lambda, pm, K, cx)
#  -.coxLoglikScore(theta, event, x, z, z2, lambda, pm, K, cx, score = FALSE)

####### logLike and score of the cox model
.coxLoglikScore = function(theta, delta, x, z, z2, lambda, pm, K, cx, score = FALSE) {
  ### theta:1:p0 for alpha, 
  ###       p0+      1 to p0+pm[1] for beta1, 
  ###       p0+pm[1]+1 to pm+pm[2] for beta2,
  ###       p0+pm[2]+1 to pm+pm[3] for beta3 etc.  
  ###
  ### x    : n x m biomarker matrix
  ### z    : covariate main effect
  ### z2   : interaction with biomarkers, 
  ###        must have an intercept term for biomarker main effect
  ### pm   : pm[i-1]+1 to pm[i] is the columns index of z2 that
  ###        interacts with biomarker i
  ###     
  m  = ncol(x)
  p0 = ncol(z) 
  p  = p0 + pm[m+1]

  beta = theta[(p0+1):p]
  gma  = theta[1:p]           ### regression coefficients
  x0   = theta[(p+1):(p+m)]   ### x0 in theta is used as cut point to Zc
  a    = cx[1, ]
  b    = cx[2, ]              
  c0   = cx[3, ]              ### c0 = cx[3, ] is used for the pLRT in loglik    
  
  Zc   = z
  bzdw = NULL
  
  for(i in 1:m) {
    idx  = (pm[i]+1):pm[i+1]
    zi   = z2[, idx]     ### use x0 as the cutpoint
    w0   = 1/(1+exp(-K*(x[, i]-x0[i])))  ### faster
    dw0  = -K*w0*(1-w0)                  ### faster
    Zc   = cbind(Zc,    zi*w0)
    bzdw = cbind(bzdw,  zi%*%beta[idx]*dw0)
  }
  zb = as.vector(Zc%*%gma)
  eb = exp(zb)
  
  S0 = cumsum(eb)             # time in decreasing order, rcumsum is not needed here
  loglik = sum(delta*(zb-log(S0))) + lambda*sum(a*log(x0/c0) + b*log((1-x0)/(1-c0)))
  if(!score) return(loglik)
  
  Zcc = cbind(Zc, bzdw)
  
  #S1  = apply(eb*Zcc, 2, cumsum) 
  S1  = matrixStats::colCumsums(eb * Zcc)
  ### This is about 23% faster using a for loop
  
  score = colSums(delta*(Zcc - S1/S0))
  
  ### add penalty term
  for(i in 1:m) score[p+i] = score[p+i]+lambda*(a[i]/x0[i]-b[i]/(1-x0[i]))
  return(list(logLik = loglik, score=score, Zc = Zc))
}

####### Cache the last theta's full loglik and score computation so that
####### fn (.coxpLoglik) and gr (.coxpLRTScore) don't duplicate work
####### when optim's L-BFGS-B calls them at the same theta (which happens
####### routinely during line search).
.coxLoglikScore_cached = local({
  last_key    = NULL
  last_result = NULL
  
  function(theta, event, x, z, z2, lambda, pm, K, cx) {
    key = theta   # extra args (event,x,z,z2,lambda,pm,K,cx) are fixed within one optim() call
    if (!is.null(last_key) && length(key) == length(last_key) &&
        isTRUE(all.equal(key, last_key, tolerance = 0))) {
      return(last_result)
    }
    res = .coxLoglikScore(theta, event, x, z, z2, lambda, pm, K, cx, score = TRUE)
    last_key    <<- key
    last_result <<- res
    res
  }
})

.coxpLRTScore = function(theta, event, x, z, z2, lambda, pm, K, cx) 
  return(-.coxLoglikScore_cached(theta, event, x, z, z2, lambda, pm, K, cx)$score)

.coxpLoglik = function(theta, event, x, z, z2, lambda, pm, K, cx)
  return(-.coxLoglikScore_cached(theta, event, x, z, z2, lambda, pm, K, cx)$logLik)

predict.coxpLRT = function(object, newdata = NULL, ...) {
  cmx = object$coxCmx
  cx0 = object$c.max
  
  if(is.null(newdata)) lp = object$linear.predictors
  else {
    #print(head(newdata))
    lp = predict(cmx, newdata)
  }
  return(lp)
}

coxpLRTcvPredErr = function(formula, formula2 = NULL, formulaList = NULL, data=list(...), 
                     c0 = 0.5, lambda = 30, K = 25, folders = 5, ...) {
  mf = model.frame(formula=formula, data=data)
  y  = model.response(mf)    # survival outcomes including censoring indicator.
  n  = length(y[, 1])
  
  index = c(0, round(seq_len(folders)*n/folders))
  J = length(index)
  tmp = rep(0, J-1)
  for (i in 1:(J-1)) {
    sel  = (index[i]+1):(index[i+1])
    dat0 = data[-sel, ]
    dat1 = data[ sel, ]
    newy =    y[ sel, ]
    
    fit  = try(coxpLRT(formula, formula2, formulaList, data = dat0, c0 = c0, lambda = lambda, K = K))
    if(is(fit, 'try-error')) {
      tmp[i] = NA
      next
    }
    cx   = fit$c.max
    cox0 = fit$coxCmx
    Z    = .getXZ(formula, formula2, formulaList, dat1, cx, K)
    Zc   = Z$Zc
    #tmp[i] = ibs(cox0, Zc, newy)
    #lp = -Zc1%*%cox0$coef
    #cdx = concordance(newy~lp)
    #tmp[i] = cdx$concordance
    theta = fit$theta
    
    cxterm = .plterm(rep(0.5, Z$m))   ### since lambda = 0, cxterm is a place holder here.
    tmp[i] = -.coxLoglikScore(theta, newy[, 2], Z$x, Z$z, Z$z2, lambda=0, Z$pm, K, cxterm, score = FALSE)
  }
  return(mean(tmp, na.rm = TRUE))
}
