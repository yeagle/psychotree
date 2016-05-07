## high-level convenience interface to mob()
wdmtree <- function(formula, data, na.action, cluster,
  control=mob_control(), ...)
{
  ## call mob
  m <- match.call(expand.dots = TRUE)
  m[[1L]] <- as.name("mob")
  m$fit <- wdmfit
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- m
  class(rval) <- c("wdmtree", class(rval))
  return(rval)
}

## methods
print.wdmtree <- function(x,
  title = "Wiener Diffusion tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

## glue code for calling wdm()
wdmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  cluster = NULL, ..., estfun = FALSE, object = FALSE, 
  wdm.alpha=NULL, wdm.tau=NULL, wdm.beta=NULL, wdm.delta=NULL)
{
  y <- RWiener::as.wiener(y)
  est <- RWiener::wdm(y, 
    alpha=wdm.alpha, tau=wdm.tau, beta=wdm.beta, delta=wdm.delta)
  res <- list(
    coefficients = est$par,
    objfun = -est$loglik,
    estfun = if(estfun) wdmscore(est) else NULL, 
    object = if(object) est else NULL
  )
  return(res)
}

# print.wienertree
# predict.wienertree
# plot.wienertree

## empirical estimation function (score function) 

CERROR <- 1e-10

wdmscore <- function(object) {
  y <- object$data[,object$yvar]

  alpha <- object$par["alpha"]
  tau <- object$par["tau"]
  beta <- object$par["beta"]
  delta <- object$par["delta"]

  n <- length(y[,1])
  res <- matrix(rep(NA,4*n), n,4)
  colnames(res) <- c("alpha", "tau", "beta", "delta")
  for (i in 1:n) {
    if (y[i,2] == "lower") {
      res[i,1] <- sclalpha(y[i,1], alpha, tau, beta, delta)
      res[i,2] <- scltau(y[i,1], alpha, tau, beta, delta)
      res[i,3] <- sclbeta(y[i,1], alpha, tau, beta, delta)
      res[i,4] <- scldelta(y[i,1], alpha, tau, beta, delta)
    }
    else if (y[i,2] == "upper") {
      res[i,1] <- sclalpha(y[i,1], alpha, tau, 1-beta, -delta)
      res[i,2] <- scltau(y[i,1], alpha, tau, 1-beta, -delta)
      res[i,3] <- sclbeta(y[i,1], alpha, tau, 1-beta, -delta)
      res[i,4] <- scldelta(y[i,1], alpha, tau, 1-beta, -delta)
    }
  }

  return(res)
}

pow <- function(x,y) x^y

kappaLT <- function(t, err=CERROR) {
  sqrt(2)*sqrt(-log(pi*err*t)/t)/pi
}

kappaST <- function(t, err=CERROR) {
  sqrt(2)*sqrt(-t*log(2*sqrt(2)*sqrt(pi)*err*sqrt(t))) + 2
}

fl01LT <- function(t, beta, kappa) {
  res <- 0
  for (k in 1:kappa) {
    res <- res +  ( k*exp(-0.5*pow(pi, 2)*pow(k, 2)*t)*sin(pi*beta*k) )
  }
  res <- pi*res
  return(res)
}

fl01ST <- function(t, beta, kappa) {
  res <- 0
  for (k in -kappa:kappa) {
    res <- res + ( (beta + 2*k)*exp(-0.5*pow(beta + 2*k, 2)/t) )
  }
  res <- 0.5*sqrt(2)*res/(sqrt(pi)*sqrt(pow(t, 3)))
  return(res)
}

fl01 <- function(t, beta) {
  kst <- kappaST(t)
  klt <- kappaLT(t)
  wlam <- kst - klt
  if(wlam < 0) fl01ST(t, beta, kst)
  else fl01LT(t, beta, klt)
}

sclalpha <- function(t, alpha, tau, beta, delta) {
  t <- t-tau

  res <- ( -beta*delta*fl01(t/alpha**2, beta) 
    * exp(-alpha*beta*delta 
    - 0.5*pow(delta, 2)*t)/pow(alpha, 2) 
    - 2*fl01(t/alpha**2, beta)*exp(-alpha*beta*delta 
    - 0.5*pow(delta, 2)*t)/pow(alpha, 3) 
    - 0 # 2*t*exp(-alpha*beta*delta - 0.5*pow(delta, 2)*t) * 0 
    / pow(alpha, 5) )

  return(res)
}

scl01tauLT <- function(t, beta, kappa) {
  res <- 0
  for (k in 1:kappa) {
    res <- res + ( 0.5*pow(pi, 2)*pow(k, 3)*exp(-0.5*pow(pi, 2)*pow(k, 2)*t)*sin(pi*beta*k) )
  }
  res <- res * pi
  return(res)
}

scl01tauST <- function(t, beta, kappa) {
  res1 <- 0
  res2 <- 0
  for (k in -kappa:kappa) {
    res1 <- ( -2*pow(beta + 2*k, 3)*exp(-pow(beta + 2*k, 2)
      / (2*t))/pow(2*t, 2) )
    res2 <- ( (beta + 2*k)*exp(-pow(beta + 2*k, 2)/(2*t)) )
  }
  res1 <- 0.5*sqrt(2)*res1/(sqrt(pi)*sqrt(pow(t, 3)))
  res2 <- 0.75*sqrt(2)*res2/(sqrt(pi)*(t)*sqrt(pow(t, 3)))
  return(res1+res2)
}

scl01tau <- function(t, beta) {
  kst <- kappaST(t)
  klt <- kappaLT(t)
  wlam <- kst - klt
  if(wlam < 0) scl01tauST(t, beta, kst)
  else scl01tauLT(t, beta, klt)
}

scltau <- function(t, alpha, tau, beta, delta) {
  t <- t-tau

  res <- ( 0.5*pow(delta, 2)*fl01(t/alpha**2, beta)
    * exp(-alpha*beta*delta - 0.5*pow(delta, 2)*t)/pow(alpha, 2) 
    - exp(-alpha*beta*delta - 0.5*pow(delta, 2)*t)
    * scl01tau((t/alpha**2), beta)
    / pow(alpha, 4) )

  return(res)
}

scl01betaLT <- function(t, beta, kappa) {
  res <- 0
  for (k in 1:kappa) {
    res <- res + ( pi*pow(k, 2)*exp(-0.5*pow(pi, 2)*pow(k, 2)*t)*cos(pi*beta*k) )
  }
  res <- res * pi
  return(res)
}

scl01betaST <- function(t, beta, kappa) {
  res <- 0
  for (k in -kappa:kappa) {
    res <- res + ( exp(-0.5*pow(beta + 2*k, 2)/t) - 0.5*(beta + 2*k)*(2*beta + 4*k)*exp(-0.5*pow(beta + 2*k, 2)/t)/t )
  }
  res <- 0.5*sqrt(2)*res/(sqrt(pi)*sqrt(pow(t, 3)))
  return(res)
}

scl01beta <- function(t, beta) {
  kst <- kappaST(t)
  klt <- kappaLT(t)
  wlam <- kst - klt
  if(wlam < 0) scl01betaST(t, beta, kst)
  else scl01betaLT(t, beta, klt)
}

sclbeta <- function(t, alpha, tau, beta, delta) {
  t <- t-tau

  res <- ( -delta*fl01(t/alpha**2, beta)*exp(-alpha*beta*delta 
    - 0.5*pow(delta, 2)*t)/alpha + exp(-alpha*beta*delta 
    - 0.5*pow(delta, 2)*t)
    * scl01beta(t, beta)
    / pow(alpha, 2) )

  return(res)
}

scldelta <- function(t, alpha, tau, beta, delta) {
  t <- t-tau

  res <- ( (-alpha*beta - delta*t)*fl01(t/alpha**2, beta)
    * exp(-alpha*beta*delta - 0.5*pow(delta, 2)*t)
    / pow(alpha, 2) )

  return(res)
}
