### Approximate function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

resboot = function(x, ...) UseMethod("resboot")

resboot.formula = function(formula, family, data=list(...), B = 100, epsilon = 0.01, ...) {
  mf = model.frame(formula=formula, data = data)
  
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
    
  control = list(B = B, epsilon = epsilon)
  if (is(y, "Surv")) {
    family = "surv";
    st = sort(y[, 1], decreasing = TRUE, index.return = TRUE)
    idx = st$ix
    y = y[idx, ]
    x = x[idx, ]
  }
      
  n.col  = ncol(x)
  if(n.col < 4) stop("Use resboot(Surv(time, event)~biomarker+trt+biomarker*trt)")
        
  # inverse cdf transformation of biomarker
  w = x[, 2]
  if((min(w) < 0) | (max(w) > 1)) {
    x[, 2] = x.cdf(x[, 2])
    transform = TRUE
  } else transform = FALSE
      
  fit = resboot.default(x, y, family, control)
  return(fit)
}

#function to calculate lik0, lik1, and LR(c) statistic
.LRatio=function(y, w, z) {
  zw  = z*w
  fit1=coxph(y~z+w+zw)
  fit0=coxph(y~z+w)
  l1=fit1$loglik[2]
  l0=fit0$loglik[2]
  LRT=2*(l1 - l0)
  return(list(lglk0 = l0, lglk1 = l1, LRT=LRT, fit1 = fit1, fit0 = fit0))
}

# max LR
.maxLR = function(y, w, z, epsilon) {
	  #epsilon = control$epsilon
	  ux = seq(0.1, 0.9, epsilon)
  
  lik = -length(w)*10
    lik0  = lik
    lglk  = lik
      lglk0 = lik
      LR    = 0
        mLR   = LR
        
        for (u in ux) {
		    w0 = ifelse(w > u, 1, 0)
	    tmp = .LRatio(y, w0, z)
	        
	        lglk0 = tmp$lglk0
	        lglk  = tmp$lglk1
		    LR    = tmp$LRT
		    
		    if(lglk > lik) {
			          lik = lglk
		          cmax = u
			        fit1 = tmp$fit1
			      }
		        
		        if(lglk0 > lik0) {
				      lik0 = lglk0
			      cmax0 = u
			            fit0 = tmp$fit0
			          }
		        
		        if(LR > mLR) mLR = LR
			  }
	  return(list(cmax0 = cmax0, cmax=cmax, mLR=mLR, fit0=fit0, fit1=fit1))
}

#bootfu return the p-value for residual bootstrap
resboot.default = function(x, y, family, control, ...) {
  # x = [intercept, biomarker, trt, interaction]
  n = nrow(x)

  w = x[, 2]
  z = x[, 3]
    
  # under \Theta0
  x0 = x[, -4]
  xs = x
    
  epsilon = control$epsilon
  B = control$B
        
  # profile Likelihood
  pflk = .maxLR(y, w, z, epsilon)
  cfit0 = pflk$fit0
  cmax0 = pflk$cmax0
	    
  cfit = pflk$fit1
  cmax = pflk$cmax
  mLRT = pflk$mLR

  x0[, 2] = ifelse(w > cmax0, 1, 0)
  xs[, 2] = ifelse(w > cmax, 1, 0)
  xs[, 4] = xs[, 2]*xs[, 3]
		  
  #Nelson Aalen cumulative baseline
  S0 = basehaz(cfit, centered = FALSE)
		  
  beta  = cfit$coeff
  beta0 = cfit0$coeff
  if (family == "surv") {
    xs = xs[, -1]
    x0 = x0[, -1]
    exb  = exp(xs%*%beta)
    exb0 = exp(x0%*%beta0)
  } 
      
  #estimated residuals under the alternative hypothesis
  chaz = S0$hazard
  tm0 = S0$time
  xchz = .appxf(chaz, tm0, y[, 1])
  uhat = (exp(-xchz))^exb
  i=1
  LRb = rep(0, B)
  while (i<=B) {
    idx = sample(1:n, replace=TRUE)
   
    uStar = uhat[idx]
    yStar = y[idx, ]
       
    tmStar=1-uStar^(1/exb0)
    yStar[, 1] = tmStar
    #print(summary(tmStar))

    pfStar = try(.maxLR(yStar, w, z, epsilon))
    if(is(pfStar, "try-error")) next
	        
    LRb[i] = pfStar$mLR
    i=i+1
  }
  pValue=sum(LRb>mLRT)/B
  fit = list(pValue = pValue, LRT.max = mLRT, c.max = cmax, cfit = cfit, LRb = LRb, B = B)
  class(fit) = 'resboot'
  return(fit)
}

print.resboot = function(x, ...) {
  pValue = x$pValue
  cat('\nResidual bootstrap p-value for interaction term = ', pValue, '\n')
  if(pValue <= 0.05) {
    cat("Cut point for biomarker using profile likelihood method = ", x$c.max, '\n')
    cat("Cox model based on this cut point:\n")
    print(x$cfit)
  }
}

plot.resboot = function(x, ...) {
  B = x$B
  max.LRT = x$LRb
  hist(max.LRT, breaks = min(20, round(B/4)))
  abline(v = x$LRT.max, lty = 2)
}
