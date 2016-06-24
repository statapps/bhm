#MCMC for Gibbs sampling
bhmGibbs<-function(x, y, family, beta, q, cx, control){
  c.n = control$c.n
  
  #generate cut point cx
  lik = thm.lik(x, y, family, beta, q, cx, control)
  if(c.n==2){
    rpt = TRUE
    D   = min(0.05, (cx[2]-cx[1])/3)
    D1  = min(0.05, cx[1]/2)
    D2  = min(0.05, (1-cx[2])/2)

    while (rpt) {
      cu = c(runif(1, cx[1]-D1, cx[1]+D), runif(1, cx[2]-D, cx[2]+D2))
      fit = thm.fit(x, y, family, cu)
      rpt = !(fit$converged)
    }
  } else {
    D = min(0.05, cx/2, (1-cx)/2)
    cu = runif(1, max(0.05, cx-D), min(cx+D, 0.95))
  }
  lik1 = thm.lik(x, y, family, beta, q, cu, control)
  alpha1 = exp(lik1 - lik)
  if( runif(1, 0, 1) < alpha1) {
    cx = cu
    lik = lik1
  }

#generate beta value using Metropolis-Hasting algorithm. 
#Candidate value is obtained from the fitted model with 
#current iteraction of cut points.
  fit = thm.fit(x, y, family, cx)
  bhat = fit$coefficient
  vb = vcov(fit)
  A = chol(vb)
  b1 = bhat + t(A)%*%rnorm(ncol(A), 0, 1)
  lik1 = thm.lik(x, y, family, b1, q, cx, control)
  alpha2 = exp(lik1 - lik)
  if(runif(1, 0, 1) < alpha2) {
    beta = b1
    lik = lik1
  }

#generate hyper parameter q
  if (c.n == 2) {
    sc1 =  log(cx[2]/(cx[2]-cx[1]))
    sc2 = -log(1-cx[2])
    q = c(1+rgamma(1,shape=1,scale=sc1), 1+rgamma(1,shape=1,scale=sc2))   
  } else {
    scale = -log(1-cx)
    q = 1+rgamma(1, 2, scale = scale)
  }
  return(list(cx=cx, beta=beta, q=q, lik=lik))
}

# fit the main threshold model
bhmFit = function(x, y, family, control){
  R   = control$R
  c.n = control$c.n
  tn  = control$thining
  var_names = colnames(x)
  
# use profile likelihood method to get initial value of cut-points
  pfit = pro.fit(x, y, family, control=list(R = 0, epsilon=0.1, c.n = c.n))

# generate initial values for parameters
  g = list(cx = pfit$c.max, beta = pfit$coefficient, q = rep(2, c.n))
  
#replication from 1 to B (burn-in)
  for (i in 1:control$B) g = bhmGibbs(x, y, family, g$beta, g$q, g$cx, control)

#replications from B+1 to R(total length of Markov Chain is B+R)
  R1 = R*tn
  bg = matrix(NaN, R, length(g$beta))
  cg = matrix(NaN, R, c.n)
  qg = matrix(NaN, R, c.n)
  i = 1
  for (j in 1:R1){
    g = bhmGibbs(x, y, family, g$beta, g$q, g$cx, control)
   
    if(j%%tn == 0){ 
      qg[i, ] = g$q
      cg[i, ] = g$cx
      bg[i, ] = g$beta
      i = i + 1
    }
  }
 
#estimates and credible interval for the thresholds
  c.max= apply(cg,2,mean)
  c.fit= thm.fit(x, y, family, c.max)
  alpha = control$alpha/2
  ptl  = c(alpha, 1-alpha)
  cqtl = apply(cg, 2, quantile, ptl)

#estimates for the regression coefficients
  coef= apply(bg,2,mean)
  coefqtl<-apply(bg,2,quantile,ptl)
   
  vcov<-cov(bg)
  if (family == "surv") var_names = var_names[-1]
  colnames(vcov) = var_names
  rownames(vcov) = var_names
  se<-sqrt(diag(vcov))

  fit = list(cg=cg,bg=bg,qg=qg,c.max=c.max,cqtl=cqtl,coefficients=coef,coefqtl=coefqtl,vcov=vcov,StdErr=se,var_names=var_names, c.fit = c.fit)
  return(fit)
}
