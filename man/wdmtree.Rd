\name{wdmtree}

\alias{wdmtree}
\alias{print.wdmtree}

\title{Wiener Diffusion Tree Models}

\description{
  Recursive partitioning based on Wiener Diffusion models.
}

\usage{
wdmtree(formula, data, na.action, cluster, control=mob_control(), ...)

%\method{print}{wdmtree}(object, \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be an object of class \code{\link[RWiener]{numdata.wiener}}
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{a data frame containing the dependent variable and
  partitioning variables.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.pass}}.}
  \item{cluster}{optional vector (typically numeric or factor) with a
    cluster ID to be employed for clustered covariances in the parameter
    stability tests.}
  \item{control}{takes \code{mob_control(\ldots)} as argument}
  \item{\dots}{Additional arguments passed to the \code{fit} function:
  \code{wdm.alpha}, \code{wdm.tau}, \code{wdm.beta}, \code{wdm.delta} (to
  fix the corresponding parameter of the Wiener Diffusion model to a
  spesific value).}
  %\item{object}{fitted model object of class \code{"wdmtree"}.}
}

\details{
  Wiener Diffusion tree models are an application of model-based recursive partitioning
  to Wiener Diffusion models. Details about the Wiener diffusion model and
  the RWiener package providing the necessary functions can be found in
  Wabersich & Vandekerckhove (2014) and the RWiener package help.
  
  Various methods are provided for \code{"wdmtree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print}).
}

\value{
  An object of S3 class \code{"wdmtree"} inheriting from class \code{"modelparty"}.
}
 
\references{
Wabersich, D., & Vandekerckhove, J. (2014). The RWiener package: An R
package providing distribution functions for the Wiener diffusion model.
The R Journal, 6(1), 49-56.
}

\examples{
## generate example dataset, fix seed for reproducability
set.seed(1)

## generate some random data
n <- 10
dat <- NULL
delta <- runif(n,-2,2)
for (i in 1:n) {
  
  if (i\%\%2 == 0) {
    age <- round(runif(1,20,30))
    delta <- 1
  }
  else {
    age <- round(runif(1,30,60))
    delta <- -1
  }

  if (i > 7) delta <- 2
  
  rw <- rwiener(100,2,.3,.5, delta)

  if(is.null(dat)) dat <- cbind(rw, age, id=as.factor(i))
  else dat <- rbind(dat, cbind(rw, age, id=as.factor(i)))
}
dat <- as.wiener(dat)
dat$group <- factor(1, levels=c(1,2))
dat$group[as.numeric(dat$id) > 7] <- factor(2, levels=c(1,2))
dat$y <- revamp(dat)

## wdmtree function call
wdmtree1 <- wdmtree(y~age+group, dat, cluster=id)

}

\keyword{tree}
