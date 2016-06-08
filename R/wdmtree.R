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

## glue code for calling wdm()
wdmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  cluster = NULL, ..., estfun = FALSE, object = FALSE, 
  wdm.alpha=NULL, wdm.tau=NULL, wdm.beta=NULL, wdm.delta=NULL)
{
  y <- RWiener::as.wiener(y)
  est <- RWiener::wdm(y, 
    alpha=wdm.alpha, tau=wdm.tau, beta=wdm.beta, delta=wdm.delta)
  res <- list(
    coefficients = est$coefficients,
    objfun = -est$loglik,
    estfun = if(estfun) estfun.wdmtree(est) else NULL, 
    object = if(object) est else NULL
  )
  return(res)
}

## estfun.wdm wrapper to aggregate scores by id, if id column is given
estfun.wdmtree <- function(x, ...) {
  res <- RWiener::estfun.wdm(x)
  if("id" %in% names(x$data)) {
    res <- cbind(res, id=x$data$id)
    res <- aggregate(. ~ id, sum, data=as.data.frame(res))[,-1]
  }
  return(res)
}

## methods
print.wdmtree <- function(x,
  title = "Wiener Diffusion tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

