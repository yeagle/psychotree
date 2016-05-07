\name{wdmtree}

\alias{wdmtree}

\alias{print.wdmtree}

\title{Wiener Diffusion Tree Models}

\description{
  Recursive partitioning based on Wiener Diffusion models.
}

\usage{
wdmtree(formula, data, na.action, cluster, \dots)

\method{print}{wdmtree}(object, \dots)
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
  \item{\dots}{Additional arguments passed to the \code{fit} function:
  \code{wdm.alpha}, \code{wdm.tau}, \code{wdm.beta}, \code{wdm.delta} (to
  fix the corresponding parameter of the Wiener Diffusion model to a
  spesific value).}
  \item{object}{fitted model object of class \code{"wdmtree"}.}
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
d1 <- rwiener(100, 2,.3,.5,2)
d2 <- rwiener(100, 2,.3,.5,0)
d3 <- rwiener(100, 2,.3,.5,-2)
## add group identification to data
d1$group <- factor('g1', levels=c('g1', 'g2', 'g3'))
d2$group <- factor('g2', levels=c('g1', 'g2', 'g3'))
d3$group <- factor('g3', levels=c('g1', 'g2', 'g3'))
## bind to one data.frame
dat <- rbind(d1,d2,d3)
## add reshaped dependent variable to dataset (y as single variable)
dat$y <- reshape.wiener(dat[,1:2])
## restructre data (columns "q" and "resp" are redundant)
dat <- dat[,c("y", "q", "resp", "group")]

## wdmtree function call
wdmtree1 <- wdmtree(y~group,dat)
# alternative call: wdmtree1 <- wdmtree(y~1|group,dat)

## wdmtree function call with alpha, tau, beta parameters fixed
wdmtree2 <- wdmtree(y~group,dat, wdm.alpha=2, wdm.tau=.3, wdm.beta=.5)
}

\keyword{tree}
