multiRoot = function(func, theta,..., verbose = FALSE, maxIter = 20, 
        thetaUp = NULL, thetaLow = NULL,
        tol = .Machine$double.eps^0.25) {
  alpha = 0.0001
  rho = 0.5
  U = func(theta, ...)
  m = length(U)
  p = length(theta)
  if (m == p) mp = TRUE  ### for m = p
  else mp = FALSE
  convergence = 0
  mU1 = sum(U^2)
  i = 1

  while(i < maxIter) {
    J = numJacobian(func, theta, m = m,...)
    if(mp) dtheta = solve(J, U)
    else {
      tJ = t(J)       ## m*1-(p*m) x (m*p) x (p*m) x (m*1)
      dtheta = solve(tJ%*%J, tJ%*%U)
    }
    theta0 = theta
    lambda = 1
    delta = 1
    Ud = sum(U*dtheta)
    ########Linear search
    while (delta > 0) {
      theta = theta0 - lambda*dtheta
      if(!is.null(thetaUp)) theta = ifelse(theta > thetaUp, thetaUp, theta)
      if(!is.null(thetaLow)) theta = ifelse(theta < thetaLow, thetaLow, theta)
      U = func(theta, ...)
      mU = sum(U^2)
      delta = mU-mU1-alpha*lambda*Ud
      lambda = rho*lambda
      #if(verbose) cat('delta = ', delta, '\n')
    }
    dU = abs(mU1 - mU)
    mU1 = mU
    i = i + 1
    if(verbose) cat("||U|| = ", mU, dU, '\n')
    if ((mU < tol) | (dU < tol)) {
      convergence = 1
      break
    }
  }
  return(list(root = theta, f.root = U, iter = i, convergence = convergence))
}
