#library(MASS)
#library(ggplot2) # used for the hazard ratio plot
#library(survival)

###  generic function for Biomarker Rectifier Models(brm) using the
###  Rectifier Linear Unit (ReLU).
brm = function(x, ...) {
  UseMethod("brm")
}
### the follow function take input like
###
###       brm(Surv(time, status) ~ w + trt + z1 + z2 + z3)  # need this
###          biomarker is the first, trt is the second, the interaction is the last
### and will fit a Cox PH model with
###
###       lambda_0(t) exp(b1*trt + b2*z1 + b3*z2 + b4*z3 +
###                        b5*max(w-c, 0) + b6*trt*max(w-c, 0))
###
brm.formula = function(formula, data=list(...), interaction = TRUE,  
                       method = c("gradient", "profile", "single"), q = 1, epsilon = NULL, ...) {
  method = match.arg(method)
  mf = model.frame(formula = formula, data = data)
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  var_names = colnames(x)

  if (class(y) == "Surv") {
    family = "surv"
    st = sort(y[, 1], decreasing = TRUE, index.return = TRUE)
    idx = st$ix
    y = y[idx, ]
    x = x[idx, ]
  }
  n.col = ncol(x)
  w = x[, 2]

  # Check for possible wrong order of formula
  if (length(unique(w)) < 4)
    stop("Either the formula is not in correct order of 'y~biomarker + trt' or the biomarker is not a continuous variable.")

  #Fit a prognostic model with biomarker term only
  if(n.col == 2) interaction = FALSE
  if (q>1) {
    cat("Since q>1, a signle indelx model will be fitted.\n")
    method = 'single'
  }
   
  # covariate name for interaction term, sort out proper variable names
  if(interaction) {
    int_names = paste(colnames(x)[2], ":", colnames(x)[3], sep="")
    var_names =  c(var_names, int_names)
    x = cbind(x, x[, 2]*x[, 3])
    colnames(x) = var_names
  }
  if(family=='surv') var_names=var_names[-1]
  #print(x[1:5, ])

  ### Use 20 points for first pass searc: seq_c for search grid points
  if(is.null(epsilon)) epsilon = (max(w) - min(w))/20;
  seq_c=seq(min(w), max(w), epsilon)

  control = list(family = family, interaction = interaction, method = method, q = q, seq_c = seq_c)
  p = length(var_names)
  if (method == 'single') {
    fit = .singleFit(x, y, control)
    for(i in 2:q) var_names[i-1] = paste('eta.', var_names[i], sep='')
    var_names[q] = 'index'
    if(interaction) var_names[p] = paste('index:', var_names[q+1], sep='')
  } else fit = brm.default(x, y, control)

  fit$call = match.call()
  fit$formula = formula
  fit$method = method
  fit$var_names = c(var_names, paste(var_names[q], '.cut', sep=""))
  return(fit)
}

brm.default = function (x, y, control, ...) {
  x = as.matrix(x)
  method = control$method
  seq_c = control$seq_c
  p = ncol(x)

  if(control$family=="surv") x=x[, -1]
  if(p == 2) x = as.matrix(x, ncol = 1)

  ### first pass search for initial value
  fit0 = .reluProFit(x, y, control)
  control$theta0 = c(fit0$coefficients, fit0$c.max)
  
  if (method == "gradient") 
    fit = .reluGradFit(x, y, control)
  else {
    lgc = as.matrix(fit0$log_c)
    
    ### second pass: refined search around possible MLE
    c2 = lgc[lgc[, 2] > (fit0$loglik - 3.84), 1]
    
    epsilon = seq_c[2]-seq_c[1]
    epsilon2 = (max(c2) - min(c2) + epsilon)/60
    control$seq_c = seq(min(c2)- 0.5*epsilon, max(c2) + 0.5*epsilon, epsilon2)
    fit = .reluProFit(x, y, control)
    tmp = rbind(lgc, fit$log_c)
    fit$log_c = tmp[order(tmp[, 1]), ]
  }
  
  ### Find s.e
  #numerical method in specified
  I_numerical = .reluInfo(fit$theta, x, y)
  I_numerical_inv = ginv(I_numerical, tol = sqrt(.Machine$double.eps))

  grdInfo = .reluGradient(fit$theta, x, y, info = TRUE)
  pi_theta = grdInfo$pi_theta
  
  #score in specified model
  sV = ginv(pi_theta) 
  
  #robust var by second derivative in misspecified
  second = grdInfo$second
  robust_second = solve(second)%*%pi_theta%*%solve(second)

  names(fit$coefficients) = colnames(x)
  fit$var = sV
  fit$linear.predictors = grdInfo$lp
  fit$first = grdInfo$U
  fit$sd.numerical = sqrt(diag(I_numerical_inv))
  fit$sd.score = sqrt(diag(sV))         
  fit$sd.robust = sqrt(diag(robust_second))     
  fit$call = match.call()
  class(fit) = "brm"
  return(fit)
}

########### Single index model #######################
.singleFit = function(x, y, control) {
  x = as.matrix(x)
  p = ncol(x)
  # fix the first index coef = 1 for model identifibility
  q = control$q-1 
  if(p<q+1) stop("Columns of x shall greater than q+1")

  w = x[, 2:(q+2)]
  w = apply(w, 2, function(x) {(x-mean(x))/sd(x)})
  # the first column of z (intercept) will be replaced by the single index
  z = x[, c(1, (q+3):p)]
  theta = rep(0, p)
  #J = .singleLik(theta, w, z, y, q)
  jmax = nlminb(theta, .singleLik, .singleGradient, w = w, z = z, y = y, q = q)
  theta = jmax$par
  beta = theta[1:(p-1)] 
  c.max = theta[p]

  a = c(1, theta[1:q])
  z[, 1] = w%*%a
  theta.s = theta[-(1:q)]
  ev = .evalLP(theta.s, z)

  #s1 = .singleGradient(theta, w, z, y, q) # check if the score function is correct.
  #s2 = .singleGradient2(theta, w, z, y, q)
  #print(rbind(s1, s2))
  info = .singleInfo(theta, w, z, y, q)
  #print(info)
  var = ginv(info)
  #print(var)
  sd.score = sqrt(diag(var))

  fit = list(coefficients = beta, var=var, sd.score=sd.score, theta = theta, 
	     iter = jmax$iterations, c.max = c.max, loglik = -jmax$objective, 
	     linear.predictors = ev$lp, y = y, index = z[, 1])
  class(fit) = 'brm'
  return(fit)
}

########### help functions ###########################
.evalLP = function(theta, x) {
  ###  Surv(time, status) ~ biomarker + trt + x2 + x3 + ...  # need this
  ### input 'x' is a nxp matrix with 2 or more columns:
  p = length(theta)
  pb= (p-1)
  w = x[, 1]    # biomarker
  
  cx = theta[p]
  ### beta = theta1, ... theta[p-1], cut point = theta[p]
  beta = theta[1:pb]
  phi = ifelse(w >= cx, w-cx, 0)
  x[, 1] = phi
  eta = ifelse(w >= cx, -beta[1], 0)
  if (length(grep(":", colnames(x)[pb])) == 1) {
    trt = x[, 2]   # trt
    x[, pb] = trt*phi
    eta = ifelse(w >= cx, -(beta[1]+beta[pb]*trt), 0)
  }
  # print(head(X))
  ### output X is a n x (bp) matrix: [phi, trt,... x_{pb-1}, phi*trt]
  lp = (x%*%beta)[, 1]
  return(list(X = x, lp = lp, eta = eta, w = w, cx = cx, phi = phi))
}

###################################################
.singleLik=function(theta, w, z, y, q) {
  a = c(1, theta[1:q])
  z[, 1] = w%*%a
  theta.s = theta[-(1:q)]
  J = .reluLglik(theta.s, z, y)
}

###################################################
.reluLglik=function(theta, x, y){
  ev = .evalLP(theta, x)
  lp  = ev$lp
  elp = exp(lp)
  
  status = y[, 2]
  ### find cost as negative of log likelihood function
  J = -sum(status*(lp-log(cumsum(elp))))
  return(J)
}

############ To be Validated for V ########################
.reluGradient = function(theta, x, y, info = FALSE){
  ev = .evalLP(theta, x)
  lp = ev$lp
  elp = exp(lp)
  status = y[, 2]
  
  Xa = cbind(ev$X, c.max = ev$eta)
  s0 = cumsum(elp)
  s1 = apply(Xa*elp, 2, cumsum)
  U = Xa - s1/s0  ### U is a n*p matrix, each row ui
  grd = -colSums(status * U)
  
  n = nrow(Xa)
  p = ncol(Xa)
  if (p>2) interaction = (length(grep(":", colnames(x)[p-1])) == 1)
  else interaction = FALSE
  # trt = X[, 2] and z2 = [1, trt]
  if (interaction) z2 = cbind(rep(1, n), ev$X[, 2]) else z2=as.matrix(rep(1, n), nrow = n)
  I = ifelse(ev$w >= ev$cx, 1, 0)
  z2_I = z2*I

  if(info) {
    mtm = function(m) {return(m%*%t(m))}
    ### code for information matrix here
    UU = apply(U, 1, mtm) ### p*p*n each with u_i*t(u_i)
    U2 = aperm(array(UU, c(p, p, n)), c(3, 1, 2))  ### now n*p*p
    pi_theta = colSums(status * U2, dims = 1)
    
    ### Find 2nd order derivative,
    ## Sm is a lower triangular matrix of 1 because time is sorted by decending
    Sm = matrix(1, n, n)
    Sm[upper.tri(Sm)] = 0
    
    XX = array(apply(Xa, 1, mtm), c(p, p, n))
    # multiply each XX(p, p, i) with elp[i], by change elp to a p*p*n array
    # with each of i th pxp matrix = elp[i]
    X2 = XX * array(rep(elp, each = p*p), c(p, p, n))
    
    ### calculate s2, a n*p*p array, X2%*%Sm like cumsum over riskset
    s2 = apply(X2, c(1, 2), function(x, y){return(y%*%x)}, Sm)
    
    ### calculate s1_i*t(s1_i)
    # aperm changes a p1*p2*p3 array to a p3*p1*p2 array with c(3, 1, 2)
    s1_2 = aperm(array(apply(s1, 1, mtm), c(p, p, n)),c(3, 1, 2))
    I_theta = -colSums(status*(s2/c(s0)-s1_2/c(s0)^2), dims = 1)
    
    #for v_2rc
    s1rc = apply(-z2_I*elp, 2, cumsum)
    V_2rc = status%*%(-z2_I - s1rc/(s0))
    
    #for v_2cc
    den=density(ev$w)
    
    pt1=which(den$x>=theta[p])[1]
    f.w.c = 0.5*(den$y[pt1]+den$y[pt1-1])

    ### beta[1] + beta[p-1]*trt 
    if(interaction) z2b = z2%*%theta[c(1, (p-1))] else z2b = z2*theta[1]
    s1cc = apply(z2b*elp, 2, cumsum)
    V_2cc = f.w.c*(status%*%(z2b-s1cc/(s0)))
    
    #file V.2
    V.2=matrix(0,nrow=p,ncol=p)
    V.2[p, p]= V_2cc
    V.2[p, 1] = V_2rc[1]
    if(interaction) V.2[p,p-1] = V_2rc[2]
    V.2[1, p] = V_2rc[1]
    if(interaction) V.2[p-1,p] = V_2rc[2]
    
    second = I_theta + V.2
    
    upi = list(U = grd, pi_theta = pi_theta, I_theta = I_theta, second = second, lp = lp)
    return(upi)
  } else {
    return(grd)
  }
}

########## main fit functions ###########
.reluProFit = function(x, y, control) {
  p = length(x[1, ])
  ### p main effect, 1 interaction, 1 threshold
  theta0 = rep(0, p+1)
  c1 = control$seq_c
  
  nc = length(c1)
  log_c = matrix(NaN, nc, 2)
  log_c[, 1] = c1
  
  logLik = -1e10
  for (j in 1:nc) {
    theta0[p+1] = c1[j]
    X = .evalLP(theta0, x)$X
    
    ### Make sure matrix X is not singular with too many no effective biomarker
    if(mean(X[, 1]>0) < 0.05) next
    fit = coxph(y~X)
    log_c[j, 2]= fit$loglik[2]
    
    if (log_c[j, 2] > logLik) {
      beta = fit$coefficients
      logLik = log_c[j, 2]
      c.max = c1[j]
    }
  }
  theta = c(beta, c.max)
  log_c = as.matrix(log_c[!is.nan(log_c[, 2]), ])
  
  if (ncol(log_c) == 1)
    log_c=t(log_c)
  else log_c = log_c
  x = x
  y = y
  time = y[, 1]
  status = y[, 2]
  
  return(list(coefficients = beta, c.max = c.max, theta = theta, loglik = logLik,
           x = x, y = y, time = time, status = status, log_c = log_c))
}

.reluGradFit = function(x, y, control) {
  theta0 = control$theta0
  lmax = nlminb(theta0, .reluLglik, .reluGradient, x = x, y = y,
                lower = c(rep(-10, 4)), upper = c(rep(10, 4)))
  theta = lmax$par
  p = length(theta)
  beta = theta[1:(p-1)]
  c.max = theta[p]
  logLik = -lmax$objective
  
  #status = y[, 2]
  ev = .evalLP(theta, x)
  return(list(coefficients = beta, c.max = c.max, theta = theta, loglik = logLik,
           y = y, index = x[, 1], iter = lmax$iterations))
}

print.brm = function(x, digits = 4, ...) {
  cat("brm: Cox PH model with ReLU\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat('Coefficients and threshold parameter\n')
  bc = cbind(unname(t(x$coefficients)), x$c.max)
  colnames(bc) = t(x$var_names) 
  print(bc, digits=digits)
}

summary.brm = function(object, alpha = 0.05,...){
  zscore = qnorm(1-alpha/2)
  sd = object$sd.score
  theta = t(cbind(unname(t(object$coefficients)), object$c.max))
  
  qtl = cbind(theta-zscore*sd, theta+zscore*sd)
  
  TAB1 = cbind(theta, sd, theta/sd, qtl[,1], qtl[,2], 2*pnorm(-abs(theta/sd)))
  colnames(TAB1) = c("Estimate", "Std.Err", "Z value", "95% CI(Low)", "95% CI(Up)", "Pr(>z)")
  rownames(TAB1) = object$var_names
  
  lp = -(object$linear.predictors)
  cidx = concordance(object$y~lp) 
  results = list(call=object$call,TAB1=TAB1, cidx = cidx)
  class(results) = "summary.brm"
  return(results)
}

print.summary.brm = function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat('Summary table of coefficients and threshold parameter\n')
  printCoefmat(x$TAB1, digits=4)
  print(x$cidx)
}

residuals.brm = function(object, type=c("Martingale", "Cox"), ...) {
  type = match.arg(type)
  ### code for martingle residuals
  x = object
  #p = length(x$x[1, ])
  #p1 = p-1
  #x1 = x$x[, 1]
  #z = x$x[, 1:p1]
  w = x$index
  
  #s0 = cumsum(exp(X%*%beta))
  status = x$y[, 2]
  elp = exp(x$linear.predictors)
  s0 = cumsum(elp)
  H0 = rev(cumsum(rev(x$status/s0))) # cumulative baseline hazard
  #H = H0*exp(X%*%beta)               # cumulative hazard
  H = H0*elp
  r = status-H                     # martingale residuals
  
  if (type=="Cox"){
    sv = survfit(Surv(H, x$status)~1)
    t=sv$time
    surv = ifelse(sv$surv==0, 0.05, sv$surv) # in case of 0
    H_e = -log(surv)
    results = list(t=t, H_e=H_e)
    class(results) = "residuals.brm"
    return(results)
  }
  
  if (type == "Martingale"){
    results = list(r=r, w=w)
    class(results) = "residuals.brm"
    return(results)
  }
}

plot.residuals.brm = function(x, type="Martingale", ...) {
  if (type=="Cox") {
    plot(x$t, x$H_e, xlab="Cox-snell residuals", ylab="Estimated cumulative hazards", 
         ylim=c(0, max(x$H_e)), xlim=c(0, max(x$t)))
    comb = cbind(x$t,x$H_e)
    comb[!is.finite(comb)] = NA
    comb = comb[complete.cases(comb),]    
    lw = lowess(comb[,2] ~ comb[,1])
    lines(lw, lwd=3, col="red")
  }
  
  if (type == "Martingale"){
    plot(x$w, x$r, xlab="Biomarker", ylab="Martingale residuals")
    lw = lowess(x$r ~ x$w)
    lines(lw, lwd=3, col="red")
  }
}

#' @export
#' @import ggplot2
#' @importFrom graphics plot
#' @title HR_plot
#' @description This function gives plot of the odds ratio.
#' @param x summary table from the fitted model.
#' @return plot odds ratio with CIs.
plot.brm = function(x, type=c("HR"), ...){
  # INPUT
  # x: TAB1=cbind(theta, sd, theta/sd, qtl[,1], qtl[,2], 2*pnorm(-abs(theta/sd))) from summary.brm
  # OUTPUT
  # the red point is the hazard ratio
  # the length of blue line is the ci
  # the number under the lower right of each line indicates which estimator
  object = summary(x)
  x = object$TAB1
  theta = x[,1]
  beta = theta[-length(theta)]
  beta_length = length(beta)
  y = exp(beta)
  
  lower = exp(x[,4][-length(theta)])
  upper = exp(x[,5][-length(theta)])
  
  df = data.frame(x=seq(1:(beta_length)), y=y, lower=lower, upper=upper)
  p = ggplot(df, aes(x = x, y=y)) + labs(x = "order of estimators", y = "odds ratio") + geom_point(colour="red")
  p= p + theme_bw()
  p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p = p + geom_hline(yintercept=1, colour = "green", lty = 3)
  p = p + geom_linerange(aes(ymin = lower, ymax = upper),colour="blue")
  p = p + geom_text(aes(label = x), hjust = 0.5, vjust = 0.5, nudge_x = -0.1, nudge_y = 0.1)
  p = p + coord_flip()
  print(p)
}

### Numerical method for information matrix
.reluInfo = function(theta, x, y){
  p = length(theta)
  h = 0.0001
  Info = matrix(0, p, p)
  for(i in 1:p) {
    theta1 = theta
    theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    Info[i, ] = .5*(.reluGradient(theta2, x, y) - .reluGradient(theta1, x, y))/h
  }
  return(Info)
}

.singleGradient = function(theta, w, z, y, q){
  a = c(1, theta[1:q])
  z[, 1] = w%*%a
  #print(head(z))
  theta.s = theta[-(1:q)]
  
  ev = .evalLP(theta.s, z)
  lp = ev$lp
  elp = exp(lp)
  status = y[, 2]

  eta = ev$eta
  w1 = w[, -1]
  Xa = cbind(w = -w1*eta, ev$X, c.max = eta)
  #print(head(Xa))
  s0 = cumsum(elp)
  s1 = apply(Xa*elp, 2, cumsum)
  U = Xa - s1/s0  ### U is a n*p matrix, each row ui
  grd = -colSums(status * U)
}

.singleGradient2 = function(theta, w, z, y, q) {
  p = length(theta) 
  h = 0.00001
  U = rep(0, p)
  for(i in 1:p) {
    theta1 = theta
    theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    U[i] = 0.5*(.singleLik(theta2, w, z, y, q) - .singleLik(theta1, w, z, y, q))/h
  }
  U
}

.singleInfo = function(theta, w, z, y, q){
  p = length(theta)
  h = 0.0001

  Info = matrix(0, p, p)
  for(i in 1:p) {
    theta1 = theta
    theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    Info[i, ] = (.singleGradient(theta2, w, z, y, q) - .singleGradient(theta1, w, z, y, q))/h 
  }
  Info = (Info+t(Info))*0.25

  return(Info)
}

###### Generate survival data ################
gendat.surv =  function(n, c0, beta, type = c("brm", "bhm")){
  type = match.arg(type)
  ### biomarkers
  x = rnorm(n, 0.2, 2)
  ### code below to generate data with ReLU
  if(type == "bhm") x1 = ifelse(x>c0, 1, 0)
  else x1 = ifelse(x >= c0, x-c0, 0)
  
  z = rbinom(n,1,0.5)  ####used binomial to determine 0 or 1####
  zx = z*x1
  X = cbind(x1, z, zx)
  h0 = .5
  h = h0*exp(X%*%beta)
  
  stime = rexp(n, h)               #Failure time.
  endstudy = runif(n, 0, 5)
  cstatus = ifelse(stime>endstudy, 0, 1) ##Censored at end of study time.
  #cat('\nCensoring: ', 1-mean(cstatus), '\n')
  stime = ifelse(stime>endstudy, endstudy, stime)
  dat = cbind(time=stime, status=cstatus, z=z, x=x)
  dat = data.frame(dat)
  return(dat)
}
