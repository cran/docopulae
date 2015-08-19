

#' Parametric Model
#'
#' \code{param} creates an initial parametric model object.
#' Unlike other model statements this function does not perform any computation.
#'
#' @param fisherIf \code{function(x, ...)}, where \code{x} is a numeric vector, usually a point from the design space.
#' It shall evaluate to the Fisher information matrix.
#' @param dDim length of \code{x}, usually the dimensionality of the design space.
#'
#' @return \code{param} returns an object of \code{class} \code{"param"}.
#' An object of class \code{"param"} is a list containing the following components:
#' \itemize{
#' \item fisherIf: argument
#' \item dDim: argument
#' \item x: a row matrix of points where \code{fisherIf} has already been evaluated.
#' \item fisherI: a list of Fisher information matrices, for each row in \code{x} respectively.
#' }
#'
#' @seealso \code{\link{fisherI}}, \code{\link{update.param}}, \code{\link{FedorovWynn}}, \code{\link{getM}}
#'
#' @example R/examples/main.R
#'
#' @export
param = function(fisherIf, dDim) {
    r = list(fisherIf=fisherIf, dDim=dDim, x=matrix(nrow=0, ncol=dDim), fisherI=list())
    class(r) = 'param'
    return(r)
}


#' Update Parametric Model
#'
#' \code{update.param} evaluates the Fisher information at uncharted points and returns an updated model object.
#'
#' @param object some model.
#' @param x a row matrix of points.
#' The number of columns shall be equal to \code{object$dDim}.
#' @param ... ignored.
#'
#' @return \code{update.param} returns an object of \code{class} \code{"param"}.
#'
#' @seealso \code{\link{param}}, \code{\link{update.desigh}}, \code{\link{update_reference}}, \code{\link{getM}}
#'
#' @examples ## see examples for param
#'
#' @export
update.param = function(object, x, ...) {
    mod = object # workaround for S3 requirement

    if (ncol(x) != mod$dDim)
        stop(paste('x shall have exactly', mod$dDim, 'columns'))

    r = mod
    x = unique(rbind(mod$x, x)) # merge x
    idcs = seq1(nrow(mod$x) + 1, nrow(x))

    if (length(idcs) != 0) {
        xx = x[idcs,, drop=F]
        r = lapply(split(xx, seq_len(nrow(xx))), mod$fisherIf)
        fisherI = c(mod$fisherI, r)
        names(fisherI) = NULL
        mod$fisherI = fisherI
    }

    mod$x = x
    return(mod)
}


## from http://adv-r.had.co.nz/Computing-on-the-language.html#substitute
substitute_q <- function(x, env) {
    call <- substitute(substitute(y, env), list(y = x))
    eval(call)
}


#' Build Density
#'
#' \code{buildf} builds the joint probabilty density given the marginal distributions and some copula.
#'
#' If \code{buildf} should build an expression, the copula shall provide distribution expressions.
#' Please note that expressions are not validated.
#'
#' @param margins either \itemize{
#' \item \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It shall return a column matrix of two, the probability densities and cumulative distributions.
#' \item list of pairs of expressions, named \code{"pdf"} and \code{"cdf"}, the probability density and cumulative distribution.
#' }
#' @param copula a copula object from package \pkg{copula}.
#' @param names (if \code{margins} is a function) a vector of names or indices, the sequence of copula parameters in \code{theta}.
#' \code{0} or \code{""} indicates copula parameters to omit.
#'
#' @return \code{buildf} returns either \itemize{
#' \item \code{function(y, theta, ...)}, the joint probability density function, if \code{margins} is a function.
#' \item the joint probabilty density as an expression, otherwise.
#' }
#'
#' @references uses substitute_q from \url{http://adv-r.had.co.nz/Computing-on-the-language.html}
#'
#' @seealso \pkg{copula}, \code{\link{expr2f}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
#'
#' @example R/examples/buildf.R
#'
#' @export
buildf = function(margins, copula, names=NULL) {
    tt = list(margins, copula, names)
    if (is.function(margins)) {
        idcs = which(names != 0 & names != '')
        if (length(idcs) == 0) {
            r = function(y, theta, ...) {
                dp = margins(y, theta, ...)
                return(copula::dCopula(dp[,2], copula) * prod(dp[,1]))
            }
        } else {
            names = names[idcs]
            r = function(y, theta, ...) {
                dp = margins(y, theta, ...)
                copula@parameters[idcs] = as.numeric(theta[names])
                return(copula::dCopula(dp[,2], copula) * prod(dp[,1]))
            }
        }
    } else { # list of expressions
        exprdist = attr(copula, 'exprdist')
        if (is.null(exprdist))
            stop('copula does not provide distribution expressions')

        margins = lapply(margins, function(margin) lapply(margin, function(e)  substitute((e), list(e=e)) )) # wrap expressions in ()
        n = paste('u', seq1(1, length(margins)), sep='')
        p = lapply(margins, function(margin) margin$cdf)
        base::names(p) = n
        d = lapply(margins, function(margin) margin$pdf)
        base::names(d) = n

        r = substitute((a)*b, list(a=substitute_q(exprdist$pdf, p), b=substitute_q(parse(text=paste(n, collapse='*'))[[1]], d)))
    }
    return(r)
}


joinLanguage = function(...) {
    x = list(...)
    x = lapply(x, function(x) if (x[[1]] == '{') as.list(x)[seq1(2, length(x))] else x)
    x = unlist(x, recursive=F, use.names=F)
    names_ = paste('x', seq1(1, length(x)), sep='')
    names(x) = names_
    return(substitute_q(parse(text=paste('{', paste(names_, collapse=';'), '}', sep=''))[[1]], x))
}

withQuotes = function(x) {
    if (is.character(x))
        return(paste('\'', x, '\'', sep=''))
    return(as.character(x))
}

#' Expression To Function
#'
#' \code{expr2f} turns an expression into \code{function(y, theta, ...)}.
#'
#' @param x an expression.
#' @param map a named list of character strings defining left assignments (\code{a="b"} => \code{a <- b}).
#' @param yMap like \code{map} with \code{a=b} resolving to \code{a <- y[b]}.
#' @param thetaMap like \code{map} with \code{a=b} resolving to \code{a <- theta[[b]]}.
#'
#' @return \code{expr2f} returns \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It evaluates expression \code{x}.
#'
#' @references uses substitute_q from \url{http://adv-r.had.co.nz/Computing-on-the-language.html}
#'
#' @seealso \code{\link{buildf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#'
#' @export
expr2f = function(x, map=NULL, yMap=NULL, thetaMap=NULL) {
    if (is.null(map))
        map = list()
    if (!is.null(yMap)) {
        yMap[] = paste('y[', lapply(yMap, withQuotes), ']', sep='')
        map = modifyList(map, yMap)
    }
    if (!is.null(thetaMap)) {
        thetaMap[] = paste('theta[[', lapply(thetaMap, withQuotes), ']]', sep='')
        map = modifyList(map, thetaMap)
    }

    map = parse(text=paste('{', paste(names(map), map, sep='=', collapse=';'), '}'))[[1]]
    return(as.function(c(alist(y=, theta=, ...=), list(joinLanguage(map, x)) )) )
}


#' Build Derivative Function for Log f
#'
#' Builds a function that evaluates to the first/second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}/\code{theta[[i]]} and \code{theta[[j]]}.
#'
#' \pkg{numDeriv} produces \code{NaN}s if the log evaluates to (negative) \code{Inf} so you may want to specify \code{logZero} and \code{logInf}.
#'
#' @param f \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' A joint probability density function.
#' @param logZero the value \code{log(f)} should return if \code{f} evaluates to \code{0}.
#' @param logInf the value \code{log(f)} should return if \code{f} evaluates to \code{Inf}.
#' @param method see \pkg{numDeriv}
#' @param method.args see \pkg{numDeriv}
#'
#' @seealso \pkg{numDeriv}, \code{\link{buildf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#'
#' @name numDerivLogf
NULL

#' @rdname numDerivLogf
#'
#' @return \code{numDerivLogf} returns \code{function(y, theta, i, ...)} which evaluates to the first derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}.
#'
#' @export
numDerivLogf = function(f, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=F)) {
    tt = list(f, logZero, logInf, method, method.args)
    f = as.function(f)

    logf = function(theta_i, y, theta, i, ...) {
        theta[i] = theta_i
        r = log(f(y, theta, ...))
        if (is.infinite(r)) # numDeriv can't handle +-Inf
            if (r == Inf)
                return(logInf)
            else
                return(logZero)
        return(r)
    }
    r = function(y, theta, i, ...) {
        return(numDeriv::grad(logf, theta[[i]], method, NULL, method.args, y, theta, i, ...))
    }
    return(r)
}


#' @rdname numDerivLogf
#'
#' @return \code{numDeriv2Logf} returns \code{function(y, theta, i, j, ...)} which evaluates to the second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]} and \code{theta[[j]]}.
#'
#' @export
numDeriv2Logf = function(f, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=F)) {
    tt = list(f, logZero, logInf, method, method.args)
    f = as.function(f)

    logf = function(theta_ij, y, theta, ij, ...) {
        theta[ij] = theta_ij
        r = log(f(y, theta, ...))
        if (is.infinite(r)) # numDeriv can't handle +-Inf
            if (r == Inf)
                return(logInf)
            else
                return(logZero)
        return(r)
    }
    if (method == 'Richardson') {
        r = function(y, theta, i, j, ...) {
            if (i == j) {
                ## sole second derivative at D[2]
                return( numDeriv::genD(logf, theta[[i]], method, method.args, y, theta, i, ...)$D[2] )
            }
            return( numDeriv::genD(logf, c(theta[[i]], theta[[j]]), method, method.args, y, theta, c(i, j), ...)$D[4] ) # sole mixed second derivative at D[4]
        }
    } else if (method == 'complex') {
        r = function(y, theta, i, j, ...) {
            r = numDeriv::hessian(logf, c(theta[[i]], theta[[j]]), method, method.args, y, theta, c(i, j), ...)
            if (i == j)
                return(r[1])
            return(r[2])
        }
    } else {
        stop('method not implemented.')
    }
    return(r)
}


#' Build Derivative Function for Log f
#'
#' Builds a function that evaluates to the first/second derivative of \code{log(f)} with respect to a predefined set of variables/variable combinations.
#'
#' While \code{numDerivLogf} relies on \pkg{numDeriv} and therefore uses finite differences to evaluate the derivatives, \code{DerivLogf} utilizes \code{Deriv} to build sub functions for each variable in \code{names}.
#' The same is true for \code{Deriv2Logf}.
#'
#' \code{Deriv} won't recognize components or parameters accessed by \code{[}, \code{[[} or \code{$} as variables (e.g. \code{theta[["beta1"]]}).
#' Therefore it's necessary to specify mappings from \code{y} and \code{theta} to the variables in \code{f}.
#'
#' @param f an expression, a joint probability density.
#' @param names a character vector of variable names.
#' @param map a named list of character strings defining left assignments (\code{a="b"} => \code{a <- b}).
#' @param yMap like \code{map} with \code{a=b} resolving to \code{a <- y[b]}.
#' @param thetaMap like \code{map} with \code{a=b} resolving to \code{a <- theta[[b]]}.
#'
#' @seealso \code{\link[Deriv]{Deriv}} in package \pkg{Deriv}, \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#' ## mind the gain regarding runtime compared to numDeriv
#'
#' @name DerivLogf
NULL

#' @rdname DerivLogf
#'
#' @return \code{DerivLogf} returns \code{function(y, theta, i, ...)} where \code{theta} is a list of parameters.
#' It evaluates to the first derivative of \code{log(f)} with respect to variable \code{i}.
#' Additionally the attribute \code{"d"} contains the list of sub functions.
#'
#' @export
DerivLogf = function(f, names, map=NULL, yMap=NULL, thetaMap=NULL) {
    logf = substitute(log((f)), list(f=f))
    d = Derivf(logf, names)
    d = lapply(d, expr2f, map, yMap, thetaMap)

    r = function(y, theta, i, ...) {
        return(d[[i]](y, theta, ...))
    }
    attr(r, 'd') = d
    return(r)
}

#' @rdname DerivLogf
#'
#' @return \code{Deriv2Logf} returns \code{function(y, theta, i, j, ...)} where \code{theta} is a list of parameters.
#' It evaluates to the second derivative of \code{log(f)} with respect to the variables \code{i} and \code{j}.
#' Additionally the attribute \code{"d2"} contains the list of sub functions.
#'
#' @export
Deriv2Logf = function(f, names, map=NULL, yMap=NULL, thetaMap=NULL) {
    logf = substitute(log((f)), list(f=f))
    d2 = Deriv2f(logf, names)
    d2 = lapply(d2, function(d2) lapply(d2, expr2f, map, yMap, thetaMap))

    r = function(y, theta, i, j, ...) {
        return(d2[[i]][[j]](y, theta, ...))
    }
    attr(r, 'd2') = d2
    return(r)
}


#' Fisher Information
#'
#' \code{fisherI} utilizes \code{nint_integrate} to evaluate the Fisher information.
#'
#' If \code{ff} is a list, it shall contain \code{dlogf} xor \code{d2logf}.
#'
#' @param ff either \itemize{
#' \item \code{function(y, theta, i, j, ...)} which evaluates to the inner part of the expectation integral/sum.
#' \item \code{list(f=function(y, theta, ...), d2logf=function(y, theta, i, j, ...))} (recommended)
#' \item \code{list(f=function(y, theta, ...), dlogf=function(y, theta, i, ...))}
#' }
#' where \code{f} is the joint probability density function and \code{dlogf}/\code{d2logf} is the first/second derivative of \code{log(f)} with respect to \code{theta[[i]]}/\code{theta[[i]]} and \code{theta[[j]]}.
#' @param theta the list of parameters.
#' @param names a vector of names or indices, the subset of parameters to use.
#' @param yspace a space, the support of \code{y}.
#' @param ... other arguments passed to \code{ff}.
#'
#' @return \code{fisherI} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{buildf}}, \code{\link{expr2f}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{nint_space}}, \code{\link{nint_transform}}, \code{\link{nint_integrate}}, \code{\link{param}}
#'
#' @examples ## see examples for param
#'
#' @export
fisherI = function(ff, theta, names, yspace, ...) {
    tt = list(ff, theta, names, yspace)

    i = 0
    j = 0
    if (inherits(ff, 'list')) {
        dlogf = ff[['dlogf']]
        d2logf = ff[['d2logf']]
        if ((is.null(dlogf) && is.null(d2logf)) || (!is.null(dlogf) && !is.null(d2logf)))
            stop('either dlogf xor d2logf shall be given')

        f = ff[['f']]
        if (!is.null(dlogf)) {
            g = function(y, ...) dlogf(y, theta, i, ...)*dlogf(y, theta, j, ...)*f(y, theta, ...)
            gd = function(y, ...) dlogf(y, theta, i, ...)**2 *f(y, theta, ...)
        } else {
            g = function(y, ...) -d2logf(y, theta, i, j, ...)*f(y, theta, ...)
            gd = g
        }
    } else {
        g = function(y, ...) ff(y, theta, i, j, ...)
        gd = g
    }

    n = length(names)
    if (n == 0)
        return(matrix(nrow=0, ncol=0))

    combs = combn(names, 2)
    r = matrix(0, nrow=n, ncol=n, dimnames=list(names, names))

    ## do off diagonal
    for (k in seq1(1, ncol(combs))) {
        i = combs[1, k]
        j = combs[2, k]

        r[i, j] = nint_integrate(g, yspace, ...)
        #print(j / ncol(combs))
        #print(rr)
    }

    ## do diagonal
    for (i in names) {
        j = i # necessary
        r[i, i] = nint_integrate(gd, yspace, ...)
        #print(j / ncol(combs))
        #print(rr)
    }

    return(mirrorMatrix(r))
}


design = function(mod, x, w, sens, args=list(), adds=list()) {
    r = list(model=mod, x=x, w=w, sens=sens, args=args, adds=adds)
    class(r) = 'desigh'
    return(r)
}

getM_ = function(m, w) {
    return(apply(sweep(m, 3, w, '*'), c(1, 2), sum))
}

sensD = function(m, Mi, ...) {
    return(apply(m, 3, function(m, Mi) sum(diag(Mi %*% m)), Mi))
}

sensDs = function(m, Mi, A, ...) {
    t1 = Mi %*% A %*% solve(t(A) %*% Mi %*% A) %*% t(A) %*% Mi
    return(apply(m, 3, function(m, t1) sum(diag(t1 %*% m)), t1))
}

getIdcs = function(names, mod) {
    if (is.numeric(names))
        return( sort(unique(names)) )

    if (length(mod$fisherI) == 0)
        stop('model shall contain at least one Fisher information matrix')
    tt = mod$fisherI[[1]]

    if (is.null(names))
        return( seq1(1, nrow(tt)) )

    nn = rownames(tt)
    if (is.null(nn)) {
        nn = colnames(tt)
        if (is.null(nn)) {
            stop('first Fisher information matrix shall contain row or column names')
        }
    }

    return(sort( unique(match(names, nn)) ))
}

getA = function(sIdcs, idcs) {
    s = length(sIdcs)
    n = length(idcs)
    i = match(sIdcs, idcs)

    r = matrix(0, nrow=n, ncol=s)
    r[(seq1(1, s) - 1)*n + i] = 1
    return(r)
}

#' Fedorov Wynn
#'
#' \code{FedorovWynn} finds a D- or Ds-optimal design using a Fedorov-Wynn-type algorithm.
#'
#' Both \code{sNames} and \code{names} default to the set of parameters for which the Fisher information is available.
#'
#' The algorithm starts from a uniform weight design.
#' In each iteration weight is redistributed to the point which has the highest sensitivity.
#' Sequence: \code{1/i}.
#' The algorithm stops when all sensitivities are below a certain absolute and relative tolerance level, or the maximum number of iterations is reached.
#'
#' @param mod some model.
#' @param sNames a vector of names or indices, the subset of parameters to optimize for.
#' @param names a vector of names or indices, the set of parameters use.
#' @param tolAbs the absolute tolerance regarding the sensitivities.
#' @param tolRel the relative tolerance regarding the sensitivities with respect to the theoretical limit.
#' @param maxIter the maximum number of iterations.
#'
#' @return \code{FedorovWynn} returns an object of \code{class} \code{"desigh"}.
#' An object of class \code{"desigh"} is a list containing the following components:
#' \itemize{
#' \item mod: argument
#' \item x: a row matrix of design points, here identical to \code{mod$x}.
#' \item w: a numeric vector of weights, for each design point respectively.
#' \item sens: a numeric vector of sensitivities, for each design point respectively.
#' \item args: a list of arguments.
#' \item adds: a list of additional (runtime) information.
#' }
#'
#' @references Fedorov, V. V. (1971) The Design of Experiments in the Multiresponse Case.
#' \emph{Theory of Probability and its Applications}, 16(2):323-332.
#'
#' Wynn, Henry P. (1970) The Sequential Generation of D-Optimum Experimental Designs.
#' \emph{The Annals of Mathematical Statistics}, 41(5):1655-1664.
#'
#' @seealso \code{\link{param}}, \code{\link{getM}}, \code{\link{reduce}}, \code{\link{plot.desigh}}, \code{\link{Defficiency}}, \code{\link{update.desigh}}
#'
#' @examples ## see examples for param
#'
#' @export
FedorovWynn = function(mod, sNames=NULL, names=NULL, tolAbs=Inf, tolRel=1e-4, maxIter=1e4) {
    args = list(FedorovWynn=list(sNames=sNames, names=names, tolAbs=tolAbs, tolRel=tolRel, maxIter=maxIter))

    if (nrow(mod$x) == 0)
        return(design(mod, mod$x, numeric(0), numeric(0), args=args))

    sIdcs = getIdcs(sNames, mod)
    idcs = getIdcs(names, mod)
    if (!all(sIdcs %in% idcs))
        stop('sNames shall be a subset of names (argument)')

    m = simplify2array(mod$fisherI)
    if (length(idcs) != dim(m)[1]) {
        m = m[idcs, idcs,]
    }

    if (anyNA(m))
        stop('Fisher information matrices shall not contain missing values')

    if ( identical(sIdcs, idcs) ) {
        sensF = sensD
        A = NULL
        target = length(idcs)
    } else {
        sensF = sensDs
        s = length(sIdcs)
        A = getA(sIdcs, idcs)
        target = s
    }

    tolAbs_ = min(tolAbs, tolRel * target)
    n = dim(m)[3]
    w = rep(1/n, n)

    for (iIter in seq1(1, maxIter)) {
        M = getM_(m, w)
        Mi = solve(M)

        sens = sensF(m, Mi, A)

        maxIdx = which(sens == max(sens))
        if (length(maxIdx) != 1)
            maxIdx = sample(maxIdx, 1)
        d = sens[maxIdx] - target
        if (d < tolAbs_)
            break

        dw = 1 / (iIter + 1)
        w = w * (1 - dw)
        w[maxIdx] = 0
        w[maxIdx] = 1 - sum(w) # equal to 'w[maxIdx] + dw'
    }

    base::names(sens) = NULL
    adds = list(FedorovWynn=list(tolBreak=d < tolAbs_, nIter=iIter))
    return(design(mod, mod$x, w, sens, args=args, adds=adds))
}


wPoint = function(x, w) {
    ## x = row matrix
    ## w = vector
    ## nrow(x) == length(w)
    return( apply(sweep(x, 1, w, '*'), 2, sum) / sum(w) )
}

#' Reduce Design
#'
#' \code{reduce} drops insignificant design points and merges design points in a certain neighbourhood.
#'
#' @param des some design.
#' @param distMax maximum euclidean distance between points to be merged.
#' @param wMin minimum weight a design point shall have to be considered significant.
#'
#' @return \code{reduce} returns an object of \code{class} \code{"desigh"}.
#' See \code{FedorovWynn} for its structural definition.
#' The sensitivity is set to \code{NA} for all design points.
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{update.desigh}}, \code{\link{plot.desigh}}
#'
#' @examples ## see examples for param
#'
#' @export
reduce = function(des, distMax, wMin=1e-6) {
    x = des$x
    w = des$w

    m = wMin <= w
    x = x[m,, drop=F]
    w = w[m]
    cl = clusterPeak(x, w, distMax)

    wg = split(w, cl)
    rx = do.call(rbind, mapply(function(x, w) wPoint(x, w), split(as.data.frame(x), cl), wg, SIMPLIFY=F))
    if (is.null(rx))
        rx = matrix(nrow=0, ncol=ncol(x))
    dimnames(rx) = NULL
    rw = vapply(wg, sum, 0)
    rw = rw / sum(rw)
    names(rw) = NULL

    ord = roworder(rx)
    rx = rx[ord,, drop=F]
    rw = rw[ord]
    rargs = des$args
    rargs$reduce = list(des=des, distMax=distMax, wMin=wMin)
    return(design(des$model, rx, rw, rep(NA, nrow(rx)), rargs, des$adds))
}


getm = function(des, mod=NULL) {
    if (is.null(mod))
        mod = des$model
    idcs = rowmatch(des$x, mod$x)
    if (anyNA(idcs))
        stop('model shall contain Fisher information matrices for each point in the design. See update_reference, update.desigh and update.model')
    return(simplify2array(mod$fisherI[idcs]))
}


#' Update Design
#'
#' \code{update.desigh} updates the underlying model and the sensitivities.
#' Usually this is necessary for custom or transformed (e.g. reduced) designs.
#'
#' @param object some design.
#' @param ... ignored.
#'
#' @return \code{update.desigh} returns an object of \code{class} \code{"desigh"}.
#' See \code{FedorovWynn} for its structural definition.
#'
#' @seealso \code{\link{reduce}}, \code{\link{update.param}}, \code{\link{Defficiency}}, \code{\link{FedorovWynn}}, \code{\link{getM}}
#'
#' @examples ## see examples for param
#'
#' @export
update.desigh = function(object, ...) {
    des = object # workaround for S3 requirement

    r = des

    mod = update(des$model, des$x)
    r$model = mod

    sIdcs = getIdcs(des$args$FedorovWynn$sNames, mod)
    idcs = getIdcs(des$args$FedorovWynn$names, mod)

    m = getm(r)
    if (length(idcs) != dim(m)[1]) {
        m = m[idcs, idcs,]
    }
    Mi = solve(getM_(m, des$w))

    if (identical(sIdcs, idcs)) {
        sens = sensD(m, Mi)
    } else {
        A = getA(sIdcs, idcs)
        sens = sensDs(m, Mi, A)
    }

    r$sens = sens
    return(r)
}


#' Update Reference Design
#'
#' \code{update_reference} updates the underlying model of some reference design.
#' The existence of Fisher information matrices for alien design points is a common requirement for comparative statistics.
#'
#' @param ref some design.
#' @param other either a single design or a list structure of designs.
#'
#' @return \code{update_reference} returns an object of \code{class} \code{"desigh"}.
#' See \code{FedorovWynn} for its structural definition.
#'
#' @seealso \code{\link{update.param}}, \code{\link{Defficiency}}
#'
#' @examples ## see examples for param
#'
#' @export
update_reference = function(ref, other) {
    other = flatten(other)
    mod = ref$model
    for (o in other) {
        mod = update(mod, o$x)
    }
    r = ref
    r$model = mod
    return(r)
}


#' Get Fisher Information
#'
#' \code{getM} returns the Fisher information corresponding to some design.
#'
#' @param des some design.
#' @param mod some model to take Fisher informations from.
#' Defaults to \code{des$model}.
#'
#' @return \code{getM} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{param}}, \code{\link{update.desigh}}, \code{\link{update.param}}
#'
#' @examples ## see examples for param
#'
#' @export
getM = function(des, mod=NULL) {
    m = getm(des, mod)
    return(getM_(m, des$w))
}


## from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
## Adds transparency to colours
add.alpha <- function(col, alpha=1){
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
}


#' Plot Design
#'
#' \code{plot.desigh} creates a one dimensional plot of sensitivities and weights.
#' If the design space has additional dimensions, the design is projected on a specified margin.
#'
#' If \code{plus=T}, \code{wDes} is specified and its sensitivities contain missing values, then the latter are linearly interpolated from the sensitivities in \code{x}.
#'
#' If \code{circles=T}, the diameter of each circle is proportional to the corresponding weight.
#'
#' @param x some design.
#' @param ... other arguments passed to plot.
#' @param margins a vector of indices, the dimensions to project on.
#' Defaults to \code{1}.
#' @param wDes a design to take weights from.
#' Defaults to \code{x}.
#' See \code{reduce}.
#' @param plus add plus symbols to the sensitivity.
#' @param circles draw weights as circles instead of as bars.
#' @param border (if drawing circles) \code{c(bottom, left, top, right)}, the relative margins to add.
#' @param sensArgs a list of arguments to use for drawing the sensitivities.
#' @param wArgs a list of arguments to use for drawing the weights.
#'
#' @references uses add.alpha from \url{http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html}
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{reduce}}
#'
#' @examples ## see examples for param
#'
#' @export
plot.desigh = function(x, ..., margins=NULL, wDes=NULL, plus=T, circles=F, border=c(0.1, 0.1, 0, 0.1), sensArgs=list(), wArgs=list()) {
    des = x # workaround for S3 requirement

    if (is.null(margins))
        margins = 1
        #margins = 1:(des$model$dDim)
    if (is.null(wDes))
        wDes = des

    if (1 < length(margins))
        stop('not yet implemented')

    ## marginal projections
    x = des$x
    sens = des$sens
    idcs = split(seq1(1, nrow(x)), lapply(margins, function(margin) x[, margin]), drop=T)
    x = x[sapply(idcs, function(idcs) idcs[1]), margins, drop=F]
    sens = sapply(idcs, function(idcs) max(sens[idcs]))

    ord = roworder(x)
    x = x[ord,, drop=F]
    sens = sens[ord]

    wx = wDes$x
    ww = wDes$w
    wsens = wDes$sens
    idcs = split(seq1(1, nrow(wx)), lapply(margins, function(margin) wx[, margin]), drop=T)
    wx = wx[sapply(idcs, function(idcs) idcs[1]), margins, drop=F]
    ww = sapply(idcs, function(idcs) sum(ww[idcs]))
    wsens = sapply(idcs, function(idcs) max(wsens[idcs]))

    args = list(...)
    if (length(margins) == 1) {
        xlim = range(x)
        if (isTRUE(circles)) {
            d = diff(xlim)
            xlim = xlim + c(-1, 1)*d*border[c(2, 4)]
        }
        if (is.null(args$ylim)) {
            ymax = max(sens)
            ylim = c(0, ymax)
            if (isTRUE(circles)) {
                d = diff(ylim)
                ylim = ylim + c(-1, 1)*d*border[c(1, 3)]
            }
        } else {
            ymax = args$ylim[2]
            ylim = NULL
        }
        xlab = colnames(x)
        if (is.null(xlab))
            xlab = paste('x[, c(', toString(margins), ')]', sep='')
        ylab = 'sensitivity'
        args = modifyList(list(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab), args)

        sensArgs = modifyList(list(col='black'), sensArgs)
        wArgs = modifyList(list(col='black'), wArgs)

        par(mar=c(5, 4, 4, 4) + 0.1)
        do.call(plot, modifyList(list(NA), args))
        do.call(lines, modifyList(list(x, sens), sensArgs))

        if (isTRUE(circles)) {
            do.call(symbols, modifyList(list(wx, rep(0, nrow(wx)), ww, inches=1/2, add=T), wArgs))
            alpha = ww/max(ww)
            idcs = which(1/256 < alpha)
            do.call(abline, modifyList(modifyList(list(v=wx[idcs], lty=3), sensArgs), list(col=add.alpha(sensArgs$col, alpha[idcs])) ))
            #points(wx, rep(0, nrow(wx)), cex=1/2)
        } else {
            if (plus) {
                alpha = ww/max(ww)
                idcs = which(1/256 < alpha)
                wx_ = wx[idcs]
                wsens_ = wsens[idcs]
                if (anyNA(wsens_)) {
                    if (nrow(x) < 2)
                        warning('need at least two design points to interpolate sensitivity')
                    else {
                        nas = is.na(wsens_)
                        wsens_[nas] = approx(x, sens, wx_[nas])$y
                    }
                }
                do.call(points, modifyList(modifyList(list(wx_, wsens_, pch='+'), sensArgs), list(col=add.alpha(sensArgs$col, alpha[idcs])) ))
            }
            par(new=T)
            do.call(plot, modifyList(list(wx, ww, type='h', xlim=args$xlim, ylim=c(0, 1), axes=F, xlab='', ylab=''), wArgs))
            do.call(axis, modifyList(list(4, col.axis=wArgs$col), wArgs))
            do.call(mtext, modifyList(list('weight', side=4, line=3), wArgs))
        }
    }
}


#' D Efficiency
#'
#' \code{Defficiency} computes the D- or Ds-efficiency measure for some design with respect to some reference design.
#'
#' Both \code{sNames} and \code{names} default to the corresponding argument which was used to find the reference design.
#'
#' D efficiency is defined as
#' \deqn{\left(\frac{\left|M(\xi,\bar{\theta})\right|}{\left|M(\xi^{*},\bar{\theta})\right|}\right)^{1/n}}{( det(M(\xi, \theta))  /  det(M(\xi*, \theta)) )**(1/n)}
#' and Ds efficiency as
#' \deqn{\left(\frac{\left|M_{11}(\xi,\bar{\theta})-M_{12}(\xi,\bar{\theta})M_{22}^{-1}(\xi,\bar{\theta})M_{12}^{T}(\xi,\bar{\theta})\right|}{\left|M_{11}(\xi^{*},\bar{\theta})-M_{12}(\xi^{*},\bar{\theta})M_{22}^{-1}(\xi^{*},\bar{\theta})M_{12}^{T}(\xi^{*},\bar{\theta})\right|}\right)^{1/s}}{( det(M11(\xi, \theta) - M12(\xi, \theta) \%*\% solve(M22(\xi, \theta)) \%*\% t(M12(\xi, \theta)))  /  det(M11(\xi*, \theta) - M12(\xi*, \theta) \%*\% solve(M22(\xi*, \theta)) \%*\% t(M12(\xi*, \theta))) )**(1/s)}
#'
#' where \eqn{M_{11}}{M11} is the submatrix corresponding to the parameters in \code{sNames}, \eqn{M_{22}}{M22} is the submatrix corresponding to the parameters in \code{names} which are not in \code{sNames}, and \eqn{M_{12}}{M12} is defined as the resulting off diagonal submatrix.
#'
#' @param des some design.
#' @param ref some design, the reference.
#' @param sNames a vector of names or indices, the subset of parameters to use.
#' @param names a vector of names or indices, the set of parameters to use.
#'
#' @return \code{Defficiency} returns a single numeric.
#'
#' @seealso \code{\link{FedorovWynn}}, \code{\link{update_reference}}, \code{\link{update.desigh}}
#'
#' @examples ## see examples for param
#'
#' @export
Defficiency = function(des, ref, sNames=NULL, names=NULL) {
    if (is.null(sNames)) {
        sNames = ref$args$FedorovWynn$sNames
    }

    if (is.null(names)) {
        names = ref$args$FedorovWynn$names
    }

    sIdcs = getIdcs(sNames, ref$model)
    idcs = getIdcs(names, ref$model)
    i = match(sIdcs, idcs)
    if (anyNA(i))
        stop('sNames shall be a subset of names (argument)')

    m = getm(des, ref$model)
    if (length(idcs) != dim(m)[1]) {
        m = m[idcs, idcs,]
    }
    M = getM_(m, des$w)

    m = getm(ref)
    if (length(idcs) != dim(m)[1]) {
        m = m[idcs, idcs,]
    }
    Mref = getM_(m, ref$w)

    if ( identical(sIdcs, idcs) ) {
        return((det(M) / det(Mref)) ** (1/length(idcs)))
    }

    s = length(sIdcs)
    M12 = M[i, -i, drop=F]
    Mref12 = Mref[i, -i, drop=F]

    return(( det(M[i, i, drop=F] - M12 %*% solve(M[-i, -i, drop=F]) %*% t(M12)) / det(Mref[i, i, drop=F] - Mref12 %*% solve(Mref[-i, -i, drop=F]) %*% t(Mref12)) ) ** (1/s))
}

