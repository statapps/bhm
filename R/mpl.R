###Fit a joint model for binary and survival clustered data with Jackknife variance
#rm(list = ls())
#library(survival)
#library(MASS)

mpl = function(formula, ...) {
  UseMethod("mpl", formula)
}

mpl.formula = function(formula, formula.glm, formula.cluster, data, weights=NULL, subset=NULL, 
                       max.iter = 300, tol = 0.005, B = 20, jackknife=FALSE, bootstrap = TRUE, ...) {
  Call = match.call()
  Call[[1]] = as.name("mpl")
  indx = match(c("formula", "formula.glm", "formula.cluster", "data", "weights", "subset", "na.action"),
               names(Call), nomatch = 0)
  
  print(Call)
  #print(indx)
  # sort the survival  data
  if (indx[1] == 0) stop("a formula argument is required")
  if(!is.null(subset)) data = subset(data, subset)

  
  mf = model.frame(formula = formula, data=data)
  s = model.response(mf)

  if(!inherits(s, "Surv")) stop("The first formula response must be a survival object")
  st = sort(s[, 1], decreasing = TRUE, index = TRUE)
  idx = st$ix
  data.sorted = data[idx, ]
    
  mf1 = model.frame(formula = formula, data=data.sorted)
  mf2 = model.frame(formula = formula.glm, data=data.sorted)
  mf3 = model.frame(formula = formula.cluster, data=data.sorted)
  cluster = mf3[, 1]
  cluster = as.factor(as.numeric(as.factor(mf3[, 1])))
  
  W.cox   = model.matrix(attr(mf1, "terms"), data=mf1)
  Z.glm   = model.matrix(attr(mf2, "terms"), data=mf2)
  
  ## remove intercept term
  W.cox = as.matrix(W.cox[, -1])

  s.cox = model.response(mf1)
  y.glm = model.response(mf2)
  
  zNames = colnames(Z.glm)
  wNames = colnames(W.cox)
  zNames = paste('glm', zNames, sep = '_')
  wNames = paste('cox', wNames, sep = '_')
  
  control = list(max.iter = max.iter, tol = tol, varsig = TRUE, B = B)
  control$varNames = c(zNames, wNames, 'sigma1', 'sigma2', 'sigma_12')
  control$weights = weights
  px = cumsum(c(length(zNames), length(wNames), 3))
  control$px = px
  control$theta0 = c(rep(0, sum(control$px[2])), 1, 1, 0)
  #control$subset = subset
  
  
  fit = mplFit(y.glm, s.cox, Z.glm, W.cox, cluster, control)
  fit$jse = NULL
  if(jackknife) {
    bootstrap = FALSE
    jfit = .mplJK(y.glm, s.cox, Z.glm, W.cox, cluster, fit$theta, control)
    fit$jse = jfit$theta.jse
  }
  if(bootstrap) {
    bfit = .mplBoot(y.glm, s.cox, Z.glm, W.cox, cluster, fit$theta, control)
    fit$jse = bfit$theta.jse
  }
  theta = fit$theta
  p = length(theta)
  fit$coefficients = theta
  fit$OR_HR = exp(theta)[1:(p-2)]
  fit$control = control
  class(fit) = 'mpl'
  return(fit)
}

.updateMPL = function(y, s, Z, W, c_matrix, U, sigma, control) {
  var_u = sigma[1]
  var_v = sigma[2]
  cov_uv= sigma[3]
  time  = s[, 1]
  event = s[, 2]
  n = length(y)

  Z1 = Z[, -1]
  #print(head(c_matrix))
  U_p = c_matrix%*%U
  ##For given U, estimate beta, gamma with survival and logistic function
  ##baseline hazard for each observation was obtained by basehaz
  logi = glm(y ~ Z1 + offset(U_p[, 1]),family=binomial)
  beta<-as.numeric(logi$coefficients) 
  se.beta<-as.numeric(summary(logi)$coefficients[,2])
  #cat('\n length s', length(s[, 1]))
  #cat('\n length W', length(W[, 1]))
  #cat('\n')
  ph = coxph(s~W+offset(U_p[, 2]))
  gamma<-as.numeric(ph$coefficients)
  se.gamma<-as.numeric(summary(ph)$coefficients[,3])
  
  tm = as.data.frame(time)
  df2<-as.data.frame(basehaz(ph, centered=FALSE))
  #print(df)
  df<-merge(tm, df2, by="time")

  lambda<-df$hazard[n:1]
  #print(lambda)

  ##For given coefficient beta and gamma
  ##update U with Newton-Ralphson (treat them as parameter)
  s_u<-vector();i_u<-vector();s_v<-vector();i_v<-vector()
  list_info<-list()
  
  flag <- 0
  for (x in 1:control$max.iter){
    cur = U
    u = U[, 1]
    v = U[, 2]
    #uv_centre<-cbind(u,v)
    U_p = c_matrix%*%U
    #u_p<-uv_p[,1];v_p<-uv_p[,2]
    
    ##calculate score and information for u
    eb = exp(Z%*%beta+U_p[, 1])
    eb1 = eb/(1+eb)
    s_u = y - eb1  
    i_u = -1*eb1*(1-eb1) 
    
    ##calculate score and information for v
    exp_all = lambda*exp(W%*%gamma+U_p[, 2])
    s_v = event - exp_all
    i_v = -1*exp_all

    ##convert the n terms of score (information) to ncentre terms and add centre-specific 
    ##score and information for bivariate normal density
    dinv = 1/(var_u*var_v-cov_uv^2)
    s_u_c<-t(c_matrix)%*%s_u-(u*var_v-v*cov_uv)*dinv; i_u_c<-t(c_matrix)%*%i_u-var_v*dinv
    s_v_c<-t(c_matrix)%*%s_v- (v*var_u-u*cov_uv)*dinv; i_v_c<-t(c_matrix)%*%i_v-var_u*dinv
    cov_12<-cov_uv/(var_u*var_v-cov_uv^2)
    A<-t(c_matrix)%*%i_u; B<-t(c_matrix)%*%i_v
    
    #uv_new = matrix(0, control$ncentre, 2)    
    ##update u and v
    for (e in 1:control$ncentre) {
      uv<-cbind(u[e],v[e])
      score<-cbind(s_u_c[e], s_v_c[e])
      info = matrix(c(i_u_c[e],cov_12,cov_12,i_v_c[e]),2,2,byrow=TRUE)
      #print(info)
      
      inverse = solve(info)
      U[e, ] = uv-score%*%inverse
      list_info[[e]] = inverse
    }
    maxu = 10
    U = ifelse(U >  maxu,  maxu, U)
    U = ifelse(U < -maxu, -maxu, U)
    esp = max(abs(cur - U))
    #print(U)
    if(esp < 0.001){ flag <- 1; break}   
  } 
  return (list(U = U, beta=beta, gamma=gamma,list_info=list_info,se.beta=se.beta,se.gamma=se.gamma,A=A,B=B))
}

####joint model utilizing all the data
mplFit = function (y, s, Z, W, centre, control) {
  #beta  = rep(0, length(Z[1, ]))
  #gamma = rep(0, length(W[1, ]))
  #sigma = c(1, 1, 0) #var_u var_v cov_uv
  #theta = c(beta, gamma, sigma)
  px = control$px
  theta = control$theta0
  beta = theta[1:px[1]]
  gamma = theta[(px[1]+1):px[2]]
  sigma = theta[(px[2]+1):px[3]]
  p     = length(theta)
  #p2    = c(p-2, p-1)
  
  ncentre = nlevels(centre)
  n = length(y)
  control$ncentre = ncentre

  var_cov<-matrix(c(sigma[1], sigma[3], sigma[3], sigma[2]),2,2)
  #print(centre)
  U = mvrnorm(ncentre, c(0,0), var_cov)

  ###create the matrix which indicates the centre ID for each patient
  c_matrix<-matrix(rep(0,n*ncentre), n, ncentre)
  for (i in 1:n){
    for (j in 1:ncentre){
      if (as.numeric(centre[i])==j) {c_matrix[i,j]=1}
    }
  }
  
  flag = 0
  for (l in 1:control$max.iter){
    theta2 = theta
    # find max penalized profile likelihood given sigma
    #pfit = .PPL(y, s, Z, W, c_matrix, U, sigma, control)
    coef1 = c(beta, gamma)
    for(ii in 1:control$max.iter) {
      coef2 = coef1
      pfit  = .updateMPL(y, s, Z, W, c_matrix, U, sigma, control)
      U = pfit$U
      coef1 = c(pfit$beta, pfit$gamma)
      espc = max(abs(coef2-coef1))
      #cat(espc, '\n')
      if(espc<control$tol) break
      #if(ii > 100) {
      # stop("No converge")}
    }
    U = pfit$U
    beta <- pfit$beta
    gamma <- pfit$gamma

    ##update variance-covariance matrix (sigma)
    l_info<-pfit$list_info
    mat_sum<-Reduce('+',l_info)
    
    matrix_uv = (t(U)%*%(U)-mat_sum)/ncentre

    #update sigma
    sigma = c(matrix_uv[1, 1], matrix_uv[2, 2], matrix_uv[1, 2])
    theta = c(beta, gamma, sigma) 
    
    # Stop iteration if difference between current and new estimates is less than tol
    if( max(abs(theta2 - theta)) < control$tol){ flag = 1; break}
  }
  if(!flag) warning("Not converge.\n") 
  se.beta  = pfit$se.beta
  se.gamma = pfit$se.gamma
  
  ####manual calcuation of var of sigma
  A = pfit$A
  B = pfit$B
  u = U[, 1]; v = U[, 2]
  a = sigma[1]; b = sigma[2]; c = sigma[3]
  if(control$varsig) {
    K<-A*B*(a*b-c^2)^2-B*(a*b-c^2)*b-A*(a*b-c^2)*a+(a*b-c^2)
    Iaa<-ncentre*b^2/(2*(a*b-c^2)^2)+1/2*sum((B*b^2+A*c^2-b)*(A*B*(2*a*b^2-2*b*c^2)-B*b^2-(2*a*b-c^2)*A+b)/(K^2))-
      1/2*sum((v^2*c^2+u^2*b^2-2*u*v*b*c)*(2*a*b^2-2*b*c^2)/((a*b-c^2)^4))
    Ibb<-ncentre*a^2/(2*(a*b-c^2)^2)+1/2*sum((B*c^2+A*a^2-a)*(A*B*(2*b*a^2-2*a*c^2)-A*a^2-(2*a*b-c^2)*B+a)/(K^2))-
      1/2*sum((u^2*c^2+v^2*a^2-2*u*v*a*c)*(2*b*a^2-2*a*c^2)/((a*b-c^2)^4))
    Icc<-ncentre*(a*b+c^2)/((a*b-c^2)^2)+sum(((B*b+A*a-1)*K-(B*b*c+A*a*c-c)*(A*B*(4*c^3-4*a*b*c)+2*B*b*c+2*A*a*c-2*c))/(K^2))-
      sum(((u^2*b+v^2*a-2*u*v*c)*(a*b-c^2)^2-(u^2*b*c+v^2*a*c-u*v*a*b-u*v*c^2)*(4*c^3-4*a*b*c))/((a*b-c^2)^4))
    Iab<-ncentre*c^2/(2*(a*b-c^2)^2)-1/2*sum(((2*B*b-1)*K-(B*b^2+A*c^2-b)*(A*B*(2*b*a^2-2*a*c^2)-A*a^2-(2*a*b-c^2)*B+a))/(K^2))+
      1/2*sum(((2*u^2*b-2*u*v*c)*(a*b-c^2)^2-(v^2*c^2+u^2*b^2-2*u*v*b*c)*(2*b*a^2-2*a*c^2))/((a*b-c^2)^4))
    Iac<-ncentre*(-1)*b*c/((a*b-c^2)^2)-1/2*sum((2*A*c*K-(B*b^2+A*c^2-b)*(A*B*(4*c^3-4*a*b*c)+2*B*b*c+2*A*a*c-2*c))/(K^2))+
      1/2*sum(((2*v^2*c-2*u*v*b)*(a*b-c^2)^2-(v^2*c^2+u^2*b^2-2*u*v*b*c)*(4*c^3-4*a*b*c))/((a*b-c^2)^4))
    Ibc<-ncentre*(-1)*a*c/((a*b-c^2)^2)-1/2*sum((2*B*c*K-(B*c^2+A*a^2-a)*(A*B*(4*c^3-4*a*b*c)+2*B*b*c+2*A*a*c-2*c))/(K^2))+
      1/2*sum(((2*u^2*c-2*u*v*a)*(a*b-c^2)^2-(u^2*c^2+v^2*a^2-2*u*v*a*c)*(4*c^3-4*a*b*c))/((a*b-c^2)^4))
    Iabc<-(-1)*matrix(c(Iaa,Iab, Iac, Iab, Ibb, Ibc, Iac, Ibc, Icc),3,3, byrow=TRUE)
    inv_Iabc<-solve(Iabc)
    # se for log(sigma_u) log(sigma_v), sigma_uv
    se.sigma = sqrt(c(inv_Iabc[1,1],inv_Iabc[2,2],inv_Iabc[3,3]))
  } 
  else { se.sigma = rep(NaN, 3) }
  theta.se = c(se.beta, se.gamma, se.sigma)
  
  ### log transoformed of sigma11 and sigma22
  #theta[p2]  = log(theta[p2])
  return(list(theta=theta, ase = theta.se))
}  

.mplJK=function (y, s, Z, W, centre, theta, control) {
  control$theta0 = theta
  n = length(centre)
  ncentre = length(unique(centre))

  control$varsig = FALSE
  
  thetai = matrix(0, ncentre, length(theta))
  w = rep(0, ncentre)
  theta.bar = 0
  p = length(theta)
  
  for (k in 1:ncentre) { 
    ####jack-knife method: each time remove one centre from the dataset
    hi = n/sum(as.numeric(centre) == k)  ### here hi = 1/wi in manuscript
    w[k] = 1/hi
    idx = (centre!=k)
    y.jk = subset(y, idx)
    s.jk = subset(s, idx)
    Z.jk = subset(Z, idx)
    W.jk = subset(W, idx)

    ###convert all centres into 1 to ncentre ID scale
    centre.jk = as.factor(as.numeric(as.factor(centre[idx])))
    
    jfit = mplFit(y.jk, s.jk, Z.jk, W.jk, centre.jk, control)
    thetai[k, ] = hi*theta - (hi-1)*jfit$theta
    theta.bar = theta.bar + w[k]*thetai[k, ]
    cat('.')
  }
  cat('\n')
  V = 0
  for (k in 1:ncentre){
    theta0 = thetai[k, ]-theta.bar
    V = V + w[k]/(1-w[k])*((theta0)%*%t(theta0))
  }
  V = V/ncentre
  theta.jse = sqrt(diag(V))
  
  return(list(theta.bar = theta.bar, theta.jse = theta.jse))
}

.mplBoot = function (y, s, Z, W, centre, theta, control) {
  n = length(centre)
  ncentre = length(unique(centre))
  control$varsig = FALSE
  
  theta.bar = 0
  B = control$B

  thetab = matrix(0, B, length(theta))
  i = 1
  while (i <= B) {
    idx = sample(n, replace = TRUE)
    y.k = y[idx]
    s.k = s[idx, ]
    Z.k = Z[idx, ]
    W.k = as.matrix(W[idx, ])
    ###convert all centres into 1 to ncentre ID scale
    centre.k = as.factor(as.numeric(as.factor(centre[idx])))
    
    st = sort(s.k[, 1], decreasing = TRUE, index = TRUE)
    ix = st$ix
    y.b = y.k[ix]
    s.b = s.k[ix, ]
    Z.b = Z.k[ix, ]
    W.b = as.matrix(W.k[ix, ])
    centre.b = centre.k[ix]

    bfit = try(mplFit(y.b, s.b, Z.b, W.b, centre.b, control))
    if (class(bfit) == "try-error") next
    thetab[i, ] = bfit$theta
    i = i + 1
    cat('.')
  }
  
  cat('\n')
  theta.jse = apply(thetab, 2, sd)
  #Vb = var(thetab)
  #theta.jse = sqrt(diag(Vb))
  theta.bar = apply(thetab, 2, mean)
  
  return(list(theta.bar = theta.bar, theta.jse = theta.jse))
}

print.mpl = function(x, digits = 3,...) {
  control = x$control
  theta = x$theta
  p = length(theta)
  
  ase = x$ase
  s.idx = (p-2):p
  cse = ase
  p2    = c(p-2, p-1)
  jse = NA*ase
  if(!is.null(x$jse)) {
    jse = x$jse
    cse = jse
  }
  
  out = cbind(Coefficients = theta, ASE = ase, JSE = jse)

  ### transfer s11, s22 to exp scale
  cse[p2] = 1/(theta[p2])*cse[p2]
  theta[p2] = log(theta[p2])

  lci = exp(theta-1.96*cse)
  uci = exp(theta+1.96*cse)
  
  out = cbind(out, OR_HR = exp(theta), Lower.CI = lci, Upper.CI = uci)

  ## NA for sigma in OR_HR column
  out[s.idx, 4] = NA
  
  #lsga = log(abs(sga))
  ### jse for log(sigma11, sigma22)
  #s.jse = jse[p2]
  
  out[p, 5] = theta[p]-1.96*cse[p]
  out[p, 6] = theta[p]+1.96*cse[p]
  pv  = 2*pnorm(-abs(theta/cse))
  pv = round(pv*10000)/10000
  out = cbind(out, pValue = pv)

  rownames(out) = control$varNames
  print(out, digits=digits)
}
