##### This is a package for the non-parametric restricted mean survival time for a 
##### continuous biomarker variable
#####

rmscb = function(x, ...) UseMethod("rmscb") 

rmscb.default=function(x, y, control, ...) {
  if (class(y) != "Surv") stop("y must be a Surv object")
  
  ### sort survival time
  #st = sort(y[, 1], decreasing = FALSE, index.return = TRUE)
  #idx = st$ix
  idx = order(y[, 1])
  y = y[idx, ]
  x = x[idx, ]

  wx = x[, 2]
  wf = ecdf(wx)
  w = wf(wx)
  
  if(is.null(control$w0)) control$w0 = seq(0.05, 0.95, 0.05)
  
  n.col = ncol(x)
  if (n.col > 3) cat("\nWarning: X has more than 3 columns, only the first two will be used.")
  ### call the oen sample function
  if (n.col == 2) fit = .rmscb1(w, y, control)
  if (n.col == 3) {
    x3 = x[, 3]
    utrt = unique(x3)
    #print(utrt)
    if(length(utrt)!=2) stop("Treatmet variable must have two levels.")
    if(min(utrt)!=1) trt = as.integer(as.factor(x3))
    fit = .rmscb2(w, y, trt, control)
  }
  #print(head(x))
  fit$w0 = control$w0
  fit$w0.original = quantile(wx, control$w0)
  class(fit) = 'rmscb'
  return(fit)
}

### Put extra parameter such as bandwidth, w0, B ... in a single list
rmsControl = function(h = 0.2, kernel = 'epan', tau = 5, B = 10, rho = 2, 
                 w0 = seq(0.1, 0.9, 0.05), sig.level = 0.95) {
  control = list(h = h, kernel = kernel, tau = tau, B = B, w0 = w0, 
                 rho = rho, sig.level=sig.level)
  return(control)
}

######### Take formula input
### fit = rmscb(Surv(time, event)~w, data = dat)
### or 
### fit = rmscb(Surv(time, event)~w+treatment, data = dat)
###
#rmscb.formula = function(formula, data=list(...), subset, tau= 5, h = 0.2, w0=NULL, 
rmscb.formula = function(formula, data, subset, na.action, tau= 5, h = 0.2, w0=NULL, 
        sig.level = 0.95, rho = 2,...){
  # mf = model.frame(formula=formula, data=data)
  cl = match.call()
  mf = match.call(expand.dots = FALSE)
  m  = match(c("formula", "data", "subset", "weights", "na.action", 
       "offset"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels = TRUE
  mf[[1L]] = quote(stats::model.frame)
  mf = eval(mf, parent.frame())
 
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  
  control = rmsControl(tau=tau, h=h, w0=w0, B=10, sig.level=sig.level, rho=rho)
  fit = rmscb.default(x, y, control)
  return(fit)
}

.rmscb1 = function(w, y, control) {
  tau = control$tau
  h   = control$h
  w0  = control$w0
  n   = length(w)
  rho = control$rho

  time = y[, 1]; status = y[, 2]
  Ind  = time>=tau

  G_fit = survfit(Surv(time, 1-status)~1)
  G = .appxf(G_fit$surv, x=G_fit$time, xout = time)
  
  if (sum(Ind)>0){
    Gtau = .appxf(G_fit$surv, x=G_fit$time, xout = tau)
    G[Ind]=Gtau
    status[Ind]=1;  time[Ind]=tau
  }
  
  ### bandwith selection if h is null
  if(is.null(h)) {
    Ind2=1:n; Ind2=Ind2[status==1]
    alpha=seq(3, 5, 0.1)
    loocv = sapply(alpha, .loocv, time=time, status=status, w=w, G=G, Idx=Ind2, rho=rho) 
    record=alpha[order(loocv)[1]]
    h=n^(-1/record)
  }                        # h = n^(-delta) ==> delta = - log(h)/log(n)
  delta = -log(h)/log(n)   #print(delta)

  ###Fit a smooth rmst function of biomarker w with IPCW method
  L = length(w0)
  m = w0; v = w0; f = w0
  for (i in 1:L) {
    kh = K_func(w, w0[i], h)
    skh = sum(kh)
    m[i]=sum(kh*time*status/G, na.rm=TRUE)/skh
    v[i]=sum(kh*(time)^2*status/G^2, na.rm=TRUE)/skh
    f[i]=skh/n
  }
  s=v-m^2

  a = control$sig.level
  ca=log(2)-log(abs(log(a)))
  
  lambda=0.6; c2=1.25; d=(2*delta*log(n))^0.5+(2*delta*log(n))^(-0.5)*(log(c2/(2*pi)))
  
  p=(n*h)^(-0.5)*(s*f*lambda)^0.5*f^(-1)*(d+ca*(2*delta*log(n))^(-0.5)) 
  
  LB = m - p; UB = m + p
  cilen = mean(2*p)
  #sum[j] = ifelse( sum(LB > bw | bw > UB) == 0, 1, 0)
  #rtmse[j]=sqrt(mean((m-bw)^2)); aloi[j]=mean(UB-LB); bandw[j]=h
  n.event = sum(status)
  fit = list(rms=m, LB=LB, UB=UB, ci.length=cilen, h=h, s=s, n=n, n.event=n.event)
}

.rmscb2 = function(w, y, trt, control) {
  tau = control$tau
  h   = control$h
  w0  = control$w0
  rho = control$rho

  time = y[, 1];        status = y[, 2]
  time1 = time[trt==1]; status1 = status[trt==1]
  time2 = time[trt==2]; status2 = status[trt==2]
  w1 = w[trt==1]; n1 = length(w1)
  w2 = w[trt==2]; n2 = length(w2); n = n1 + n2 

  G_fit1 = survfit(Surv(time1, 1-status1)~1)
  G_fit2 = survfit(Surv(time2, 1-status2)~1)
  
  G1 = .appxf(G_fit1$surv, x=G_fit1$time, xout = time1)
  G2 = .appxf(G_fit2$surv, x=G_fit2$time, xout = time2)
  Ind1 = time1>=tau
  Ind2 = time2>=tau
  
  if(sum(Ind1)>0){
    Gtau1 = .appxf(G_fit1$surv, x=G_fit1$time, xout = tau)
    G1[Ind1]=Gtau1
    status1[Ind1]=1; time1[Ind1]=tau
  }
  
  if(sum(Ind2)>0){
    Gtau2 = .appxf(G_fit2$surv, x=G_fit2$time, xout = tau)
    G2[Ind2]=Gtau2
    status2[Ind2]=1; time2[Ind2]=tau
  }

  ### bandwith selection if h is null
  if(is.null(h)) {
    alpha=seq(3, 5, 0.1)
    Idx1=1:n1; Idx1=Idx1[status1==1]
    Idx2=1:n2; Idx2=Idx2[status2==1]

    loocv1 = sapply(alpha, .loocv, time=time1, status=status1, w=w1, G=G1, Idx=Idx1, rho=rho)
    loocv2 = sapply(alpha, .loocv, time=time2, status=status2, w=w2, G=G2, Idx=Idx2, rho=rho)
    
    loocv = loocv1 + loocv2
    record=alpha[order(loocv)[1]]
    h=n^(-1/record)
  }                        # h = n^(-delta) ==> delta = - log(h)/log(n)
  delta = -log(h)/log(n)   #print(delta)

  ### fit the smooth curves
  m1 = w0; m2 = w0; v1 = w0; v2 = w0; f = w0
  for (i in 1:length(w0)) {
    kh=K_func(w1, w0[i], h) 
    m1[i]=sum(kh*time1*status1/G1, na.rm=TRUE)/sum(kh)
    v1[i]=sum(kh*time1^2*status1/G1^2, na.rm=TRUE)/sum(kh)
    
    kh=K_func(w2, w0[i], h) 
    m2[i]=sum(kh*time2*status2/G2, na.rm=TRUE)/sum(kh)
    v2[i]=sum(kh*time2^2*status2/G2^2, na.rm=TRUE)/sum(kh)
    
    kh=K_func(w, w0[i], h) 
    f[i]=mean(kh)
  }
  m = m1-m2
  s1= v1-m1^2; s2=v2-m2^2
  s = s1*n/n1+s2*n/n2

  a = control$sig.level; ca=log(2)-log(abs(log(a)))
  lambda=0.6; c2=1.25; d=(2*delta*log(n))^0.5+(2*delta*log(n))^(-0.5)*(log(c2/(2*pi)))
  p=(n*h)^(-0.5)*(s*f*lambda)^0.5*f^(-1)*(d+ca*(2*delta*log(n))^(-0.5)) 
  
  LB = m - p; UB = m + p
  cilen = mean(2*p)
  #sum[j] = ifelse( sum(LB > bw | bw > UB) == 0, 1, 0)
  #rtmse[j]=sqrt(mean((m-bw)^2)); aloi[j]=mean(UB-LB); bandw[j]=h
  n.event = sum(status)
  fit = list(rms=m, LB=LB, UB=UB, ci.length=cilen, h=h, s=s, n=n, n.event=n.event)
}

plot.rmscb = function(x, x2=NULL, xlab = 'Biomarker', ylab='RMST',...){
  rms = x$rms; w0 = x$w0; LB = x$LB; UB = x$UB
  ymin = min(LB) - .25
  ymax = max(UB) + .25

  df = cbind(w0, rms, LB, UB)
  if(!is.null(x2)) {
    rms2 = x2$rms; LB2 = x2$LB; UB2 = x2$UB
    df = cbind(w0, rms, LB, UB, rms2, LB2, UB2)
  }
  df = data.frame(df)
  p = ggplot(df, aes(w0, rms)) + geom_line() + labs(x=xlab, y=ylab)
  p = p + geom_ribbon(aes(ymin=LB, ymax=UB, alpha=0.2), show.legend=FALSE, fill='yellow')
  if(!is.null(x2)) {
    p = p + geom_line(aes(y=rms2))
    p = p + geom_ribbon(aes(ymin=LB2, ymax=UB2, alpha=0.2), show.legend=FALSE)
  }
  #p = p + geom_ribbon(aes(ymin=LB, ymax=UB, fill='gray70'))

  #plot(w0, rms, ylim = c(ymin, ymax), type = 'n', xlab = 'biomarker')
  #lines(w0, rms, lwd = 3, col = 'black')
  
  # scb
  #lines(w0, LB, lwd = 2, lty = 2, col = 'black')
  #lines(w0, UB, lwd = 2, lty = 2, col = 'black')
  #abline(h=0, lwd = 3)
  print(p)
  p
}

print.rmscb = function(x, ...){
  out = summary(x)
  print(out)
}

summary.rmscb = function(object, ...) {
  rms = object$rms; w0 = object$w0; LB = object$LB; UB = object$UB
  results = data.frame(cbind(w0, rms, LB, UB))
}

### Kernel function
K_func<-function(w, w0, h, kernel = c("epanechnikov", "gaussian", "rectangular", "triangular", "biweight", "cosine", "optcosine")) {
  kernel = match.arg(kernel)
  x = w-w0
  ax = abs(x)
  esp = 1e-40

  kh = switch(kernel, gaussian = ifelse(ax < 5*h, dnorm(x, sd = h), esp), 
         rectangular = ifelse(ax < h, 0.5/h, esp), 
         triangular = ifelse(ax < h, (1 - ax/h)/h, esp),
         epanechnikov = ifelse(ax < h, 3/4 * (1 - (ax/h)^2)/h, esp), 
         biweight = ifelse(ax < h, 15/16 * (1 - (ax/h)^2)^2/h, esp),
         cosine = ifelse(ax < h, (1 + cos(pi * x/h))/(2*h), esp),
         optcosine = ifelse(ax < h, pi/4 * cos(pi * x/(2*h))/h, esp)
         )
  return(kh)
}

### Leave one out cross validation
.loocv = function(a, time, status, w, G, Idx, rho) {
  n  = length(w)
  ms = Idx
  h = n^(-1/a)
  for (k in 1:length(Idx)) {
    j=Idx[k]
    kh=K_func(w[-j], w[j], h)  ### K_func(w, w0, h, kernel) default is "epanechnikov"
    mh=sum(kh*time[-j]*status[-j]/G[-j], na.rm=TRUE)/sum(kh)
    ms[k]=(abs(time[j]-mh))^rho
  }
  return(sum(ms))
}

### x-fold cross validation 
.xfcv = function(a, time, status, w, G, rho, fold = 5) {
  n  = length(w)
  size = n/fold
  index = seq(0, n, size)
  J = length(index)
  if (index[J] != n) cat("Warning: n/K is not an integer.\n")
    
  h = n^(-1/a)
  mse = 0
  for(i in 1:(J-1)) {
    sel = (index[i]+1):(index[i+1])
      
    ### Select training set, X[-sel, ] means "sel" rows will be removed.
    time.train = time[-sel]
    status.train = status[-sel]
    w.train = w[-sel]
    G.train = G[-sel]
      
    ### select testing set
    time.test = time[sel]
    status.test = status[sel]
    w.test = w[sel]
      
    ### Fit model and find beta, the code can be replaced when 
    ### a different model is needed.
    #fit.train = lm(y.train~X.train)
    for(j in 1:size) {
      if(status.test[j] == 1) {
        kh = K_func(w.train, w.test[j], h)  ### K_func(w, w0, h, kernel) default is "epanechnikov"
        mh = sum(kh*time.train*status.train/G.train, na.rm=TRUE)/sum(kh)
        mse = mse + (abs(time.test[j] - mh))^rho
      }
    }
  }
  return(sum(mse))
}

getdat = function(n, censoring_pt, shape=10) {
  w = runif(n,0,1)
  fail_time = rweibull(n, shape=shape, scale = w)
  
  #censoring = runif(n,0,censoring_pt)
  censoring = rexp(n, rate = censoring_pt)
  status = ifelse(fail_time<=censoring, 1, 0)
  cat('censoring = ', (1-mean(status))*100, '%\n') 
  
  time = ifelse(fail_time <= censoring, fail_time, censoring)

  #idx  = order(time)
  #w    = w[idx]
  #status = status[idx]
  #time   = time[idx]
  
  data = cbind(time, status, w)
  return(data.frame(data))
}

getdat2 = function(n, censoring_pt) {
  w = runif(n,0,1)
  fail_time = rweibull(n, shape=0.5, scale = w)
  
  #censoring = runif(n,0,censoring_pt)
  censoring = rexp(n, rate = censoring_pt)
  status = ifelse(fail_time<=censoring, 1, 0)
  cat('censoring = ', (1-mean(status))*100, '%\n') 
  
  time = ifelse(fail_time <= censoring, fail_time, censoring)
  data = cbind(time, status, w)
  return(data.frame(data))
}

