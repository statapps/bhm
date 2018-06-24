bhm = function(x, ...) UseMethod("bhm")

bhm.default = function(x, y, family, control, ...) {
  cat("bhm: Biomarker thershold models\n")
  x = as.matrix(x)
  method = control$method

  if(method == 'Bayes')
    fit = bhmFit(x, y, family, control)
  else
    fit = prolikFit(x, y, family, control)

  fit$family = family
  fit$control = control
  fit$call = match.call()
  class(fit) = "bhm"
  return(fit)
}

bhm.formula = function(formula, family, data=list(...), control=list(...), ...){
  mf = model.frame(formula=formula, data=data)
  
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  
  if (class(y) == "Surv") {
    family = "surv";
    st = sort(y[, 1], decreasing = TRUE, index.return = TRUE)
    idx = st$ix
    y = y[idx, ]
    x = x[idx, ]
  }
  
  control = do.call("bhmControl", control)
  n.col  = ncol(x)

  # inverse cdf transformation of biomarker
  w = x[, 2]
  # Check for possible wrong order of formula
  if (length(unique(w)) < 4)
    stop("Either the formula is not in correct order of 'y~biomarker + trt' or the biomarker is not a continuous variable.")

  if((min(w) < 0) | (max(w) > 1)) {
    x[, 2] = x.cdf(x[, 2])
    transform = TRUE
  } else transform = FALSE


  #Fit a prognostic model with biomarker term only
  if(n.col == 2) {
    control$interaction = FALSE
  }

  # covariate name for interaction term
  if(control$interaction) {
    int_names = paste(colnames(x)[2], ":", colnames(x)[3], sep="")

    # Check if there is a main effect for biomarker variable
    if(control$biomarker.main){
      var_names =  c(colnames(x)[1:3], int_names)
      x1 = cbind(x[, 1:3], x[, 2]*x[, 3])
      if(n.col>3) {
        var_names = c(var_names, colnames(x)[4:n.col])
        x1 = cbind(x1, x[, 4:n.col])
      }
      x = x1
      colnames(x) = var_names
    } else {
      colnames(x)[2] = int_names
    }
  }
  #print(x[1:5, ])

  n.col = length(x[1, ])
  if(family == "surv"){
    n.col = n.col - 1
  }
  
  if (length(control$sigma0.inv) == 1 & n.col > 1)
    control$sigma0.inv = diag(rep(control$sigma0.inv, n.col))

  if (length(control$beta0) == 1)
    control$beta0 = rep(control$beta0, n.col)

  fit = bhm.default(x, y, family, control, ...)

  # transoform the value of biomarker back to original value
  if(transform) {
    fit$c.max = quantile(w, fit$c.max)
    qtlName = rownames(fit$cqtl)
    fit$cqtl = matrix(quantile(w, fit$cqtl), nrow = 2)
    rownames(fit$cqtl) = qtlName
  }

  fit$call = match.call()
  fit$formula = formula
  return(fit)
}

bhmControl=function(method = 'Bayes', interaction = TRUE, biomarker.main = TRUE, alpha = 0.05, B=50, R=100, thin = 1, epsilon = 0.01, c.n = 1, beta0=0, sigma0 = 10000) {

  if(method != 'profile' && method != 'Bayes')
    stop("Please use either 'Bayes' or 'profile' method for model fit")
  if (!is.numeric(B) || B <= 0) 
    stop("value of 'B' must be > 0")
  if ((!is.numeric(R) || R <= 0)&&method == 'Bayes')
    stop("number of replication 'R' must be > 0")

  if (!is.logical(biomarker.main))
    stop("'biomarker.main' must be a logical value: TURE/FALSE") 
  if (!is.logical(interaction))
    stop("'interaction' must be a logical value: TURE/FALSE")
  if (interaction == FALSE && biomarker.main == FALSE) {
    cat('Interaction = FALSE, biomarker.main effect reset to TRUE\n')
    biomarker.main = TRUE
  }

  if (!is.numeric(alpha) || alpha <= 0||alpha>=1) 
    stop("number of replication 'alpha' must be between 0 and 1")
  if (c.n > 2 || c.n < 1)
    stop("number of cutpoints 'c.n' must be either 1 or 2")
  if (!is.numeric(sigma0) || sigma0 <= 0) 
    stop("value of 'sigma' [varince for beta prior] must be > 0")
  sigma0.inv = solve(sigma0)

  return(list(method = method,  interaction = interaction, biomarker.main = biomarker.main, B=B, R=R, thin = thin, epsilon = epsilon, alpha = alpha, c.n=c.n, beta0 = beta0, sigma0.inv = sigma0.inv))
}

summary.bhm=function(object,...){
  x = object
  family = x$family
  c.max = x$c.max
  
  if (family == "binomial") {
    TAB1<-t(rbind(Estimate=x$coefficients,StdErr=x$StdErr,CredibleInterval=x$coefqtl))
    rownames(TAB1) = x$var_names
    TAB2<-t(rbind(Estimate=x$c.max,CredibleInterval=x$cqtl))

    #threshold parameter
    if (length(TAB2[, 1]) == 2) {
      rownames(TAB2)=c("lower","higer")
    } else {
      rownames(TAB2)=c("cut point")
    }
  }
  
  if (family == "surv") {
    TAB1<-t(rbind(Estimate=x$coefficients,StdErr=x$StdErr,CredibleInterval=x$coefqtl))
    TAB2<-t(rbind(Estimate=x$c.max,CredibleInterval=x$cqtl))
  }
  results = list(call=x$call,TAB1=TAB1, TAB2=TAB2, c.fit=x$c.fit, c.max = x$c.max, var_names = x$var_names)

  class(results)<-"summary.bhm"
  return(results)
}

print.summary.bhm<-function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nRegression coefficients:\n")
 
 
  printCoefmat(x$TAB1,digits=4,P.values=FALSE)

  bname=x$var_names[1]
  cat('\n', bname, 'biomarker threshold:\n')
  printCoefmat(x$TAB2,digits=4,P.values=FALSE)
  c.max = round(x$c.max*10000)/10000
  cat('\nConditional regression coefficients given', bname, 'biomarker = ', c.max, '\n')
  print(summary(x$c.fit))
}

print.bhm = function(x,...) {
   c.max = round(x$c.max*10000)/10000
   bname = x$var_names[1]
   cat("Call:\n")
   print(x$call)
   cat("\nCoefficients:\n")
   print(x$coefficients)
   cat("\n", bname, "Thresholds:\n")
   print(c.max)
   cat('\nConditional regression coefficient given ', bname, 'biomarker = ', c.max, '\n')  
   print(x$c.fit)
}

plot.bhm = function(x, type = c("profile", "density"), ...) {
  type = match.arg(type)
  cgrid = x$cgrid
  lgrid = x$lgrid
  if(type == "profile") {
    if(x$control$c.n == 1) {
      plot(cgrid, lgrid, type = 'n', xlab = 'Biomarker', ylab = 'Log likelihood')
      lines(cgrid, lgrid, lwd = 2)
      title("Profile log likelihood")
    }
    else cat("Profile plot for c.n = 1 only. \n")
  }
  if(type == "density") {
    plot(x$cg)
  }
}
