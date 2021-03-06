% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nint.R
\name{nint_transform}
\alias{nint_transform}
\title{Transform Integral}
\usage{
nint_transform(f, space, trans, funcDimToF = 0, zeroInf = 0)
}
\arguments{
\item{f}{\code{function(x, ...)}, an integrand.}

\item{space}{a space or list structure of spaces.}

\item{trans}{a list of named lists, each containing \code{dIdcs}, \code{g} and \code{giDgi} or \code{giDg}, where\itemize{
\item \code{dIdcs} is an integer vector of indices, the dimensions to transform
\item \code{g=function(x[dIdcs])} mapping \code{x[dIdcs]} to \code{y}
\item \code{giDgi=function(y)} returning a list of two, the inverse \code{gi(y) = x[dIdcs]} and the first derivatives of \code{gi(y)} with respect to \code{y}
\item or \code{giDg=function(y)} returning the inverse and the first derivatives of \code{g(x[dIdcs])} with respect to \code{x[dIdcs]}.
}}

\item{funcDimToF}{an integer vector of indices, the dimensions to look for function dimensions to transform to interval dimensions.
\code{0} indicates all dimensions.}

\item{zeroInf}{a single value, used when \code{f} returns \code{0} and the Jacobian is infinite.}
}
\value{
\code{nint_transform} returns either a named list containing the transformed integrand and space, or a list of such.
}
\description{
\code{nint_transform} applies monotonic transformations to an integrand and a space or list structure of spaces.
Common use cases include the probability integral transform, the transformation of infinite limits to finite ones and function dimensions to interval dimensions.
}
\details{
Interval dimensions and function dimensions returning interval dimensions only.

If a transformation is vector valued, that is \code{y = c(y1, ..., yn) = g(c(x1, ..., xn))}, then each component of \code{y} shall exclusively depend on the corresponding component of \code{x}.
So \code{y[i] = g[i](x[i])} for an implicit function \code{g[i]}.

The transformation of function dimensions to interval dimensions is performed after the transformations defined by \code{trans}.
Consecutive linear transformations, \code{g(x[dIdx]) = (x[dIdx] - d(x)[1])/(d(x)[2] - d(x)[1])} where \code{d} is the function dimension at dimension \code{dIdx}, are used.
Deciding against this transformation probably leads to considerable loss in computational performance.
}
\examples{
library(mvtnorm)
library(SparseGrid)

dfltNCube = nint_integrateNCube


## 1D, normal pdf
mu = 137
sigma = mu/6
f = function(x) dnorm(x, mean=mu, sd=sigma)
space = nint_space(nint_intvDim(-Inf, Inf))

tt = nint_transform(f, space,
                    list(nint_tanTransform(mu + 3, sigma*1.01, dIdcs=1)))
tt$space
ff = Vectorize(tt$f); curve(ff(x), tt$space[[1]][1], tt$space[[1]][2])

nint_integrate(tt$f, tt$space) # returns 1


## 2D, normal pdf

## prepare for SparseGrid
ncube = function(dimension)
    SparseGrid::createIntegrationGrid('GQU', dimension, 7) # rather sparse!
ncube = nint_integrateNCube_SparseGrid(ncube)
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', ncube, envir=environment(nint_integrate))

mu = c(1, 2)
sigma = matrix(c(1, 0.7,
                 0.7, 2), nrow=2)
f = function(x) {
    if (all(is.infinite(x))) # dmvnorm returns NaN in this case
        return(0)
    return(dmvnorm(x, mean=mu, sigma=sigma))
}

# plot f
x1 = seq(-1, 3, length.out=51); x2 = seq(-1, 5, length.out=51)
y = outer(x1, x2, function(x1, x2) apply(cbind(x1, x2), 1, f))
contour(x1, x2, y, xlab='x[1]', ylab='x[2]', main='f')

space = nint_space(nint_intvDim(-Inf, Inf),
                   nint_intvDim(-Inf, Inf))

tt = nint_transform(f, space,
                    list(nint_tanTransform(mu, diag(sigma), dIdcs=1:2)))
tt$space

# plot tt$f
x1 = seq(tt$space[[1]][1], tt$space[[1]][2], length.out=51)
x2 = seq(tt$space[[2]][1], tt$space[[2]][2], length.out=51)
y = outer(x1, x2, function(x1, x2) apply(cbind(x1, x2), 1, tt$f))
contour(x1, x2, y, xlab='x[1]', ylab='x[2]', main='tt$f')

nint_integrate(tt$f, tt$space) # doesn't return 1
# tan transform is inaccurate here

# probability integral transform
dsigma = diag(sigma)
t1 = list(g=function(x) pnorm(x, mean=mu, sd=dsigma),
          giDg=function(y) {
              x = qnorm(y, mean=mu, sd=dsigma)
              list(x, dnorm(x, mean=mu, sd=dsigma))
          },
          dIdcs=1:2)

tt = nint_transform(f, space, list(t1))

# plot tt$f
x1 = seq(tt$space[[1]][1], tt$space[[1]][2], length.out=51)
x2 = seq(tt$space[[2]][1], tt$space[[2]][2], length.out=51)
y = outer(x1, x2, function(x1, x2) apply(cbind(x1, x2), 1, tt$f))
contour(x1, x2, y, xlab='x[1]', ylab='x[2]', main='tt$f')

nint_integrate(tt$f, tt$space) # returns almost 1


## 2D, half sphere
f = function(x) sqrt(1 - x[1]^2 - x[2]^2)
space = nint_space(nint_intvDim(-1, 1),
                   nint_funcDim(function(x)
                        nint_intvDim(c(-1, 1)*sqrt(1 - x[1]^2))))

# plot f
x = seq(-1, 1, length.out=51)
y = outer(x, x, function(x1, x2) apply(cbind(x1, x2), 1, f))
persp(x, x, y, theta=45, phi=45, xlab='x[1]', ylab='x[2]', zlab='f')

tt = nint_transform(f, space, list())
tt$space

# plot tt$f
x1 = seq(tt$space[[1]][1], tt$space[[1]][2], length.out=51)
x2 = seq(tt$space[[2]][1], tt$space[[2]][2], length.out=51)
y = outer(x1, x2, function(x1, x2) apply(cbind(x1, x2), 1, tt$f))
persp(x1, x2, y, theta=45, phi=45, xlab='x[1]', ylab='x[2]', zlab='tt$f')

nint_integrate(tt$f, tt$space) # returns almost 4/3*pi / 2


## 2D, constrained normal pdf
f = function(x) prod(dnorm(x, 0, 1))
space = nint_space(nint_intvDim(-Inf, Inf),
                   nint_funcDim(function(x) nint_intvDim(-Inf, x[1]^2)))

tt = nint_transform(f, space, list(nint_tanTransform(0, 1, dIdcs=1:2)))

# plot tt$f
x1 = seq(tt$space[[1]][1], tt$space[[1]][2], length.out=51)
x2 = seq(tt$space[[2]][1], tt$space[[2]][2], length.out=51)
y = outer(x1, x2, function(x1, x2) apply(cbind(x1, x2), 1, tt$f))
persp(x1, x2, y, theta=45, phi=45, xlab='x[1]', ylab='x[2]', zlab='tt$f')

nint_integrate(tt$f, tt$space) # Mathematica returns 0.716315


assign('nint_integrateNCube', dfltNCube, envir=environment(nint_integrate))
}
\seealso{
\code{\link{nint_integrate}}, \code{\link{nint_space}}, \code{\link{nint_tanTransform}}, \code{\link{fisherI}}
}
