

#' Parametric Model
#'
#' \code{param} creates an initial parametric model object.
#' Unlike other model statements this function does not perform any computation.
#'
#' @param fisherIf \code{function(x, ...)}, where \code{x} is a vector, usually a point from the design space.
#' It shall evaluate to the Fisher information matrix.
#' @param dDim length of \code{x}, usually the dimensionality of the design space.
#'
#' @return \code{param} returns an object of \code{class} \code{"param"}.
#' An object of class \code{"param"} is a list containing at least the following components:
#' \itemize{
#' \item fisherIf: argument
#' \item x: a row matrix of points where \code{fisherIf} has already been evaluated.
#' \item fisherI: a list of Fisher information matrices, for each row in \code{x} respectively.
#' }
#'
#' @seealso \code{\link{fisherI}}, \code{\link{update.param}}, \code{\link{Dsensitivity}}, \code{\link{getM}}, \code{\link{Defficiency}}
#'
#' @example R/examples/main.R
#'
#' @export
param = function(fisherIf, dDim) {
    r = list(fisherIf=fisherIf, x=matrix(0, nrow=0, ncol=dDim), fisherI=list())
    class(r) = 'param'
    return(r)
}


#' Update Parametric Model
#'
#' \code{update.param} evaluates the Fisher information at uncharted points and returns an updated model object.
#'
#' @param object a model.
#' @param x either a row matrix of points or a design, or a list structure of matrices or designs.
#' The number of columns/the dimensionality of the design space shall be equal to \code{ncol(object$x)}.
#' @param ... ignored.
#'
#' @return \code{update.param} returns an object of \code{class} \code{"param"}.
#' See \code{\link{param}} for its structural definition.
#'
#' @seealso \code{\link{param}}, \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
update.param = function(object, x, ...) {
    mod = object # workaround for S3 requirement

    x = flatten(x)
    x = lapply(x, function(x) if (inherits(x, 'desigh')) x$x else x)
    x = do.call(rbind, x)

    if (ncol(x) != ncol(mod$x))
        stop(paste('x shall have exactly', ncol(mod$x), 'columns'))

    x = unique(rbind(mod$x, x)) # merge x
    idcs = seq1(nrow(mod$x) + 1, nrow(x))
    r = mod

    if (length(idcs) != 0) {
        xx = x[idcs,, drop=F]
        fisherI = lapply(split(xx, seq_len(nrow(xx))), mod$fisherIf)
        fisherI = c(mod$fisherI, fisherI)
        names(fisherI) = NULL
        r$fisherI = fisherI
        r$x = x
    }

    return(r)
}


#' Build probability density or mass Function
#'
#' \code{buildf} builds a joint probability density or mass function from marginal distributions and a copula.
#'
#' Please note that expressions are not validated.
#'
#' If \code{continuous} is \code{FALSE}, dimensionality shall be 2 and both dimensions shall be discrete.
#' The joint probability mass is defined by
#' \deqn{C(F_{1}(y_{1}),F_{2}(y_{2}))-C(F_{1}(y_{1}-1),F_{2}(y_{2}))-C(F_{1}(y_{1}),F_{2}(y_{2}-1))+C(F_{1}(y_{1}-1),F_{2}(y_{2}-1))}{C(F1(y[1]), F2(y[2])) - C(F1(y[1] - 1), F2(y[2])) - C(F1(y[1]), F2(y[2] - 1)) + C(F1(y[1] - 1), F2(y[2] - 1))}
#' where \eqn{C}, \eqn{F_{1}}{F1}, and \eqn{F_{2}}{F2} depend on \eqn{\theta} and \eqn{y_{i}\ge0}{y[i] \ge 0}.
#'
#' @param margins either \itemize{
#' \item \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' It shall return a column matrix of two, the probability densities and cumulative distributions.
#' \item a list of pairs of expressions, each named \code{"pdf"} and \code{"cdf"}, the probability density and cumulative distribution.
#' }
#' @param continuous \code{TRUE} if margins are continuous. See details.
#' @param copula if \code{margins} is \itemize{
#' \item a function then either a copula object from package \pkg{copula} or \code{function(u, theta, ...)}, a probability density function if \code{continuous} else a cumulative distribution function.
#' \item a list then either a copula object from package \pkg{copula} which contains distribution expressions or an expression for the probability density if \code{continuous} else the cumulative distribution which uses \code{u1},\code{u2},...
#' }
#' @param parNames if \itemize{
#' \item (optional) \code{margins} is a function and \code{copula} is a copula object then a vector of names or indices, the sequence of copula parameters in \code{theta}.
#' \code{0} or \code{""} identifies copula parameters to skip.
#' \item \code{margins} is a list and \code{copula} is a copula object then a named list of names or indices, mapping parameters in \code{theta} to copula parameter variables.
#' See \code{copula@exprdist}.
#' }
#' @param simplifyAndCache (if \code{margins} is a list) simplify and cache the result using \code{\link[Deriv]{Simplify}} and \code{\link[Deriv]{Cache}} from package \pkg{Deriv} if available.
#'
#' @return \code{buildf} returns \code{function(y, theta, ...)}, the joint probability density or mass function.
#'
#' @seealso \pkg{copula}, \code{\link[Deriv]{Simplify}}, \code{\link[Deriv]{Cache}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
#'
#' @example R/examples/buildf.R
#'
#' @export
buildf = function(margins, continuous, copula, parNames=NULL, simplifyAndCache=T) {
    tt = list(margins, copula, parNames)

    if (is.list(margins)) { # symbolic
        if (inherits(copula, 'copula')) { # copula object
            att = attr(copula, 'exprdist')
            if (!is.null(att)) {
                expr = if (continuous) att$pdf else att$cdf

                if (length(copula@parameters) != 0) {
                    if (is.null(parNames))
                        stop('parNames shall be specified in this use case')
                    gets = lapply(parNames, function(n) substitute(theta[[n]], list(n=n)))
                    names(gets) = names(parNames)
                    expr = substituteDirect(expr, gets)
                }

                copula = expr
            } # else: next stop
        }

        if (!is.language(copula) && !identical(copula, 1)) # special case for indepCopula
            stop('copula shall be a copula object containing distribution expressions or an expression for the probability density in this use case')

        ui = paste('u', seq1(1, length(margins)), sep='')

        cdfs = lapply(margins, function(margin) in.brackets(margin[['cdf']]))
        names(cdfs) = ui

        if (continuous) {
            pdfs = lapply(margins, function(margin) in.brackets(margin[['pdf']]))
            names(pdfs) = ui

            uProd = parse(text=paste(ui, collapse='*'))[[1]]

            f = substitute((a)*b, list(a=substituteDirect(copula, cdfs),
                                       b=substituteDirect(uProd, pdfs)))

        } else {
            cdfs1 = lapply(margins, function(margin)
                in.brackets(substituteDirect(margin[['cdf']],list(y=quote((y + c(-1, 0)))))))
            names(cdfs1) = ui

            cdfs2 = lapply(margins, function(margin)
                in.brackets(substituteDirect(margin[['cdf']], list(y=quote((y + c(0, -1)))))))
            names(cdfs2) = ui

            cdfs12 = lapply(margins, function(margin)
                in.brackets(substituteDirect(margin[['cdf']], list(y=quote((y + c(-1, -1)))))))
            names(cdfs12) = ui

            c0 = substituteDirect(copula, cdfs)
            c1 = substitute(if (y[1] == 0) 0 else { cc }, list(cc=substituteDirect(copula, cdfs1)))
            c2 = substitute(if (y[2] == 0) 0 else { cc }, list(cc=substituteDirect(copula, cdfs2)))
            c12 = substitute(if (any(y == 0)) 0 else { cc }, list(cc=substituteDirect(copula, cdfs12)))

            f = substitute((c0)-(c1)-(c2)+(c12), list(c0=c0, c1=c1, c2=c2, c12=c12))
        }

        if (simplifyAndCache) {
            if (requireNamespace('Deriv', quietly=T)) {
                f = Deriv::Cache(Deriv::Simplify(Deriv::deCache(f)))
            }
        }

        r = as.function(c(alist(y=, theta=, ...=), list(f)))
        return(r)
    }

    if (!is.function(margins))
        stop('margins shall either be a list or a function')

    if (inherits(copula, 'copula')) {
        C = copula
        CopulaF = if (continuous) copula::dCopula else copula::pCopula
        Idcs = which(parNames != 0 & parNames != '')

        if (length(Idcs) == 0) {
            copula = function(u, theta, ...) CopulaF(u, C, ...)

        } else {
            ParNames = parNames[Idcs]
            copula = function(u, theta, ...) {
                C@parameters[Idcs] = as.numeric(theta[ParNames])
                return(CopulaF(u, C, ...))
            }
        }
    }

    if (!is.function(copula))
        stop('copula shall be a copula object or a probability density function in this use case')

    if (continuous) {
        r = function(y, theta, ...) {
            dp = margins(y, theta, ...)
            return(copula(dp[,2], theta) * prod(dp[,1]))
        }

    } else {
        r = function(y, theta, ...) {
            dp = margins(y, theta, ...)
            dp1 = margins(y + c(-1, 0), theta, ...)
            dp2 = margins(y + c(0, -1), theta, ...)
            dp12 = margins(y + c(-1, -1), theta, ...)
            pp = rbind(dp[,2], dp1[,2], dp2[,2], dp12[,2])
            r = copula(pp, theta) %*% c(1, -1, -1, 1)
            return(r)
        }
    }

    return(r)
}


in.brackets = function(x) {
    r = substitute((a), list(a=x))
    return(r)
}


#' Build Derivative Function for Log f
#'
#' \code{numDerivLogf}/\code{numDeriv2Logf} builds a function that evaluates to the first/second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}/\code{theta[[i]]} and \code{theta[[j]]}.
#'
#' \pkg{numDeriv} produces \code{NaN}s if the log evaluates to (negative) \code{Inf} so you may want to specify \code{logZero} and \code{logInf}.
#'
#' @param f \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' A joint probability density function.
#' @param isLogf set to \code{TRUE} if \code{f} is already \code{log(f)}.
#' @param logZero the value \code{log(f)} should return if \code{f} evaluates to \code{0}.
#' @param logInf the value \code{log(f)} should return if \code{f} evaluates to \code{Inf}.
#' @param method,side,method.args see \code{\link[numDeriv]{grad}} and \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}.
#'
#' @seealso \code{\link[numDeriv]{grad}} and \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}, \code{\link{buildf}}, \code{\link{DerivLogf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#'
#' @name numDerivLogf
NULL

#' @rdname numDerivLogf
#'
#' @details \code{numDerivLogf} passes \code{method}, \code{side} and \code{method.args} directly to \code{numDeriv::grad}.
#'
#' @return \code{numDerivLogf} returns \code{function(y, theta, i, ...)} which evaluates to the first derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}.
#'
#' @export
numDerivLogf = function(f, isLogf=FALSE, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', side=NULL, method.args=list()) {
    tt = list(f, logZero, logInf, method, side, method.args)
    f = as.function(f)
    if (isLogf)
        log = function(x) x

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
        return(numDeriv::grad(logf, theta[[i]], method, side, method.args, y, theta, i, ...))
    }
    return(r)
}


#' @rdname numDerivLogf
#'
#' @details \code{numDeriv2Logf} duplicates the internals of \code{numDeriv::hessian} to gain speed.
#' The defaults for \code{method.args} are \code{list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)}.
#'
#' @return \code{numDeriv2Logf} returns \code{function(y, theta, i, j, ...)} which evaluates to the second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]} and \code{theta[[j]]}.
#'
#' @export
numDeriv2Logf = function(f, isLogf=FALSE, logZero=.Machine$double.xmin, logInf=.Machine$double.xmax/2, method='Richardson', method.args=list()) {
    tt = list(f, logZero, logInf, method, method.args)
    f = as.function(f)
    if (isLogf)
        log = function(x) x
    method.args_ = list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
    method.args_[names(method.args)] = method.args

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
            if (i == j)
                return( numDeriv::genD(logf, theta[[i]], method, method.args_, y, theta, i, ...)[['D']][2] ) # sole second derivative at D[2]
            return( numDeriv::genD(logf, c(theta[[i]], theta[[j]]), method, method.args_, y, theta, c(i, j), ...)[['D']][4] ) # sole mixed second derivative at D[4]
        }
    } else if (method == 'complex') {
        r = function(y, theta, i, j, ...) {
            r = numDeriv::hessian(logf, c(theta[[i]], theta[[j]]), method, method.args_, y, theta, c(i, j), ...)
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
#' \code{DerivLogf}/\code{Deriv2Logf} builds a function that evaluates to the first/second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}/\code{theta[[i]]} and \code{theta[[j]]}.
#'
#' While \code{numDerivLogf} relies on the package \pkg{numDeriv} and therefore uses finite differences to evaluate the derivatives, \code{DerivLogf} utilizes the package \pkg{Deriv} to build sub functions for each parameter in \code{parNames}.
#' The same is true for \code{Deriv2Logf}.
#'
#' @param f \code{function(y, theta, ...)}, where \code{theta} is a list of parameters.
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' @param preSimplify simplify the body of \code{f} using functions from package \pkg{Deriv}.
#' @param ... other arguments passed to \code{\link[Deriv]{Deriv}} from package \pkg{Deriv}.
#'
#' @seealso \pkg{Deriv}, \code{\link[Deriv]{Deriv}} in package \pkg{Deriv}, \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{fisherI}}
#'
#' @examples ## see examples for param
#' ## mind the gain regarding runtime compared to numDeriv
#'
#' @name DerivLogf
NULL

prepareDeriv = function(f, parNames, preSimplify) {
    f = body(f)
    logf = substitute(log(f), list(f=f))
    if (preSimplify)
        logf = Deriv::Cache(Deriv::Simplify(Deriv::deCache(logf)))

    x = parNames
    names(x) = rep('theta', length(x))

    return(list(logf=logf, x=x))
}

#' @rdname DerivLogf
#'
#' @return \code{DerivLogf} returns \code{function(y, theta, i, ...)} which evaluates to the first derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]}.
#' The attribute \code{"d"} contains the list of sub functions.
#'
#' @export
DerivLogf = function(f, parNames, preSimplify=T, ...) {
    prep = prepareDeriv(f, parNames, preSimplify)

    expr = Derivf(prep$logf, prep$x, ...)
    d = lapply(expr, function(expr)
        as.function(c(alist(y=, theta=, ...=), list(expr))))
    names(d) = names(expr)

    r = function(y, theta, i, ...) {
        return(d[[i]](y, theta, ...))
    }
    attr(r, 'd') = d
    return(r)
}

#' @rdname DerivLogf
#'
#' @return \code{Deriv2Logf} returns \code{function(y, theta, i, j, ...)} which evaluates to the second derivative of \code{log(f(y, theta, ...))} with respect to \code{theta[[i]]} and \code{theta[[j]]}.
#' The attribute \code{"d2"} contains the list of sub functions.
#'
#' @export
Deriv2Logf = function(f, parNames, preSimplify=T, ...) {
    prep = prepareDeriv(f, parNames, preSimplify)

    expr = Deriv2f(prep$logf, prep$x, ...)
    d2 = lapply(expr, function(expr) {
        r = lapply(expr, function(expr)
            as.function(c(alist(y=, theta=, ...=), list(expr))))
        names(r) = names(expr)
        return(r)
    })
    names(d2) = names(expr)

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
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' @param yspace a space, the support of \code{y}.
#' @param ... other arguments passed to \code{ff}.
#'
#' @return \code{fisherI} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{buildf}}, \code{\link{numDerivLogf}}, \code{\link{DerivLogf}}, \code{\link{nint_space}}, \code{\link{nint_transform}}, \code{\link{nint_integrate}}, \code{\link{param}}
#'
#' @examples ## see examples for param
#'
#' @export
fisherI = function(ff, theta, parNames, yspace, ...) {
    tt = list(ff, theta, parNames, yspace)

    i = 0
    j = 0
    if (inherits(ff, 'list')) {
        dlogf = ff[['dlogf']]
        d2logf = ff[['d2logf']]
        if ((is.null(dlogf) && is.null(d2logf)) || (!is.null(dlogf) && !is.null(d2logf)))
            stop('either dlogf xor d2logf shall be specified')

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

    n = length(parNames)
    if (n == 0)
        return(matrix(nrow=0, ncol=0))

    combs = combn(parNames, 2)
    r = matrix(0, nrow=n, ncol=n, dimnames=list(parNames, parNames))

    ## do off diagonal
    for (k in seq1(1, ncol(combs))) {
        i = combs[1, k]
        j = combs[2, k]

        r[i, j] = nint_integrate(g, yspace, ...)
    }

    ## do diagonal
    for (i in parNames) {
        j = i # necessary
        r[i, i] = nint_integrate(gd, yspace, ...)
    }

    return(mirrorMatrix(r))
}


#' Design
#'
#' \code{design} creates a custom design object.
#'
#' @param x a row matrix of points.
#' @param w a vector of weights.
#' Length shall be equal to the number of rows in \code{x} and sum shall be equal to \code{1}.
#' @param tag a list containing additional information about the design.
#'
#' @return \code{design} returns an object of \code{class} \code{"desigh"}.
#' An object of class \code{"desigh"} is a list containing at least this function's arguments.
#'
#' @seealso \code{\link{Wynn}}, \code{\link{reduce}}, \code{\link{getM}}, \code{\link{plot.desigh}}, \code{\link{Defficiency}}, \code{\link{update.param}}
#'
#' @examples ## see examples for param
#'
#' @export
design = function(x, w, tag=list()) {
    if (length(w) != nrow(x))
        stop('length of w shall be equal to the number of rows in x')
    if (!isTRUE(all.equal(sum(w), 1))) # R forced me to write this
        stop('weights shall sum to 1')

    r = list(x=x, w=w, tag=tag)
    class(r) = 'desigh'
    return(r)
}


getDAPar = function(mod, A, parNames) {
    # A = NULL || chr || int || row named mat || mat
    # parNames = NULL || chr || int

    # gather fI, fNames
    fI = NULL
    fNames = NULL

    if (length(mod$fisherI) != 0) {
        fI = mod$fisherI[[1]]

        fNames = rownames(fI)
        if (is.null(fNames)) {
            fNames = colnames(fI)
        }
    }

    # parNames
    if (is.character(parNames)) {
        if (!is.null(fNames)) {
            parNames = match(parNames, fNames)

            if (anyNA(parNames))
                stop('parNames doesn\'t correspond to the names in the first Fisher information matrix')
        }
    }
    # = NULL || chr || int

    # A
    if (is.null(A) || (is.matrix(A) && is.null(rownames(A)))) # D, raw D_A
        return(list(parNames=parNames, A=A))

    n = NULL
    aIdcs = NULL
    if (is.character(A)) {
        n = A
    } else if (is.vector(A)) {
        aIdcs = A
    } else {
        n = rownames(A)
    }

    if (!is.null(n)) { # resolve parNames
        if (is.character(parNames)) { # parNames = chr
            n2 = parNames
        } else { # parNames = NULL || int
            if (is.null(fNames))
                stop('in this case the first Fisher information matrix shall have names')

            if (is.null(parNames)) { # parNames = NULL
                n2 = fNames
            } else { # parNames = int
                n2 = fNames[parNames]
            }
        }

        aIdcs = match(n, n2)
        if (anyNA(aIdcs))
            stop('parameter names in A don\'t correspond to the subset of parameters')
    }

    if (is.null(parNames)) { # parNames = NULL
        if (is.null(fI))
            stop('in this case the model shall contain at least one Fisher information matrix')

        k = nrow(fI)
    } else { # parNames = chr || int
        k = length(parNames)
    }

    if (is.character(A) || is.vector(A)) { # D_s
        s = length(A)
        A = matrix(0, nrow=k, ncol=s)
        A[(seq1(1, s) - 1)*k + aIdcs] = 1
        return(list(parNames=parNames, A=A))
    }

    # named D_A
    nA = matrix(0, nrow=k, ncol=ncol(A))
    nA[aIdcs,] = A
    return(list(parNames=parNames, A=nA))
}

getm = function(mod, x, parNames=NULL) {
    if (nrow(x) == 0)
        return(array(dim=c(length(parNames), length(parNames), 0)))

    xidcs = rowmatch(x, mod$x)
    if (anyNA(xidcs))
        stop('model shall contain Fisher information matrices for each point. See update.param')

    fi = mod$fisherI[[xidcs[1]]]
    if (is.null(parNames) ||
        isTRUE(all.equal(parNames, rownames(fi))) || isTRUE(all.equal(parNames, colnames(fi))) ||
        isTRUE(all.equal(parNames, seq1(1, nrow(fi)))))
        return(simplify2array(mod$fisherI[xidcs]))

    return(vapply(xidcs, function(xidx) mod$fisherI[[xidx]][parNames, parNames], matrix(0, nrow=length(parNames), ncol=length(parNames))))
}

getM_ = function(m, w) {
    # fastest solution by now
    d = dim(m)
    nrow = d[1]
    ncol = d[2]
    r = matrix(0.0, nrow=nrow, ncol=ncol)
    jj = 1:ncol
    for (i in 1:nrow) {
        for (j in jj)
            r[i, j] = sum(m[i, j,]*w)
    }
    return(r)
}


Dsensitivity_anyNAm = 'Fisher information matrices shall not contain missing values'

#' D Sensitivity
#'
#' \code{Dsensitivity} builds a sensitivity function for the D-, D_s or D_A-optimality criterion which relies on defaults to speed up evaluation.
#' \code{Wynn} for instance requires this behaviour/protocol.
#'
#' Indices and rows of an unnamed matrix supplied to argument \code{A} correspond to the subset of parameters defined by argument \code{parNames}.
#'
#' For efficiency reasons the returned function won't complain about \emph{missing arguments} immediately, leading to strange errors.
#' Please ensure that all arguments are specified at all times.
#' This behaviour might change in future releases.
#'
#' @param A for \itemize{
#' \item D-optimality: \code{NULL}
#' \item D_s-optimality: a vector of names or indices, the subset of parameters of interest.
#' \item D_A-optimality: either \itemize{
#'   \item directly: a matrix without row names.
#'   \item indirectly: a matrix with row names corresponding to the parameters.
#'   }
#' }
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' Defaults to the parameters for which the Fisher information is available.
#' @param defaults a named list of default values.
#' The value \code{NULL} is equivalent to absence.
#'
#' @return \code{Dsensitivity} returns \code{function(x=NULL, desw=NULL, desx=NULL, mod=NULL)}, the sensitivity function.
#' It's attributes contain this function's arguments.
#'
#' @references E. Perrone & W.G. Müller (2016) Optimal designs for copula models, Statistics, 50:4, 917-929, DOI: 10.1080/02331888.2015.1111892
#'
#' @seealso \code{\link{docopulae}}, \code{\link{param}}, \code{\link{Wynn}}, \code{\link{plot.desigh}}
#'
#' @examples ## see examples for param
#'
#' @encoding UTF-8
#' @export
Dsensitivity = function(A=NULL, parNames=NULL, defaults=list(x=NULL, desw=NULL, desx=NULL, mod=NULL)) {
    # tr(M^-1 A (A^T M^-1 A)^-1 A^T M^-1 m)  = M, A, m
    # 3 M = sum(desw*desm)
    # 2   desm = (mod, desx, parNames)
    # 1     parNames = (A, parNames)
    # 1 A = (A, parNames)
    # 2 m = (mod, x, parNames)
    # 1   parNames = (A, parNames)
    dx = defaults[['x']]
    ddesw = defaults[['desw']]
    ddesx = defaults[['desx']]
    dmod = defaults[['mod']]

    # 1
    d1 = F
    if (!is.null(dmod)) {
        tt = try(getDAPar(dmod, A, parNames), silent=T)
        if (class(tt) != 'try-error') {
            d1 = T
            dparNames = tt$parNames
            dA = tt$A
        } else
            warning(tt)
    }

    # 2
    d2_1 = F
    d2_2 = F
    if (d1) {
        if (!is.null(ddesx)) {
            d2_1 = T
            ddesm = getm(dmod, ddesx, dparNames)
            if (anyNA(ddesm))
                stop(Dsensitivity_anyNAm)
        }

        if (!is.null(dx)) {
            d2_2 = T
            dm = getm(dmod, dx, dparNames)
            if (anyNA(dm))
                stop(Dsensitivity_anyNAm)
            dm = matrix(dm, nrow=prod(dim(dm)[1:2])) # prepare I
        }
    }

    # 3
    d3 = F
    if (d2_1 && !is.null(ddesw)) {
        d3 = T
        dM = getM_(ddesm, ddesw)
        dMi = solve(dM)
        if (is.null(dA))
            dt1 = dMi
        else
            dt1 = dMi %*% dA %*% solve(t(dA) %*% dMi %*% dA) %*% t(dA) %*% dMi
        dt1 = c(t(dt1)) # prepare I
    }

    # 4
    if (d3 && d2_2)
        dr = c(dt1 %*% dm) # I: tr(dt1 %*% dm[,,i])

    r = function(x=NULL, desw=NULL, desx=NULL, mod=NULL) {
        # 1
        u1 = F
        if (is.null(mod))
            mod = dmod
        else
            u1 = T

        if (u1) {
            u1 = T
            tt = getDAPar(mod, A, parNames)
            parNames = tt$parNames
            A = tt$A
        } else {
            parNames = dparNames
            A = dA
        }

        # 2
        u2_1 = F
        if (is.null(desx))
            desx = ddesx
        else
            u2_1 = T

        if (u1 || u2_1) {
            u2_1 = T
            desm = getm(mod, desx, parNames)
            if (anyNA(desm))
                stop(Dsensitivity_anyNAm)
        } else
            desm = ddesm

        u2_2 = F
        if (is.null(x))
            x = dx
        else
            u2_2 = T

        if (u1 || u2_2) {
            u2_2 = T
            m = getm(mod, x, parNames)
            if (anyNA(m))
                stop(Dsensitivity_anyNAm)
            m = matrix(m, nrow=prod(dim(m)[1:2])) # prepare I
        } else
            m = dm

        # 3
        u3 = F
        if (is.null(desw))
            desw = ddesw
        else
            u3 = T

        if (u2_1 || u3) {
            u3 = T
            M = getM_(desm, desw)
            Mi = solve(M)
            if (is.null(A))
                t1 = Mi
            else
                t1 = Mi %*% A %*% solve(t(A) %*% Mi %*% A) %*% t(A) %*% Mi
            t1 = c(t(t1)) # prepare I
        } else
            t1 = dt1

        # 4
        if (u3 || u2_2)
            r = c(t1 %*% m) # I: tr(dt1 %*% dm[,,i])
        else
            r = dr

        return(r)
    }
    attributes(r) = list(A=A, parNames=parNames, defaults=defaults)

    return(r)
}


#' Wynn
#'
#' \code{Wynn} finds an optimal design using a sensitivity function and a Wynn-algorithm.
#'
#' See \code{\link{Dsensitivity}} and it's return value for a reference implementation of a function complying with the requirements for \code{sensF}.
#'
#' The algorithm starts from a uniform weight design.
#' In each iteration weight is redistributed to the point which has the highest sensitivity.
#' Sequence: \code{1/i}.
#' The algorithm stops when all sensitivities are below a specified tolerance level or the maximum number of iterations is reached.
#'
#' @param sensF \code{function(x=NULL, desw=NULL, desx=NULL, mod=NULL)}, a sensitivity function.
#' It's attribute \code{"defaults"} shall contain identical \code{x} and \code{desx}, and \code{sensF(desw=w)} shall return sensitivities corresponding to each point in \code{x}.
#' @param tol the tolerance level regarding the sensitivities.
#' @param maxIter the maximum number of iterations.
#'
#' @return \code{Wynn} returns an object of \code{class} \code{"desigh"}.
#' See \code{\link{design}} for its structural definition.
#'
#' @references Wynn, Henry P. (1970) The Sequential Generation of D-Optimum Experimental Designs.
#' \emph{The Annals of Mathematical Statistics}, 41(5):1655-1664.
#'
#' @seealso \code{\link{Dsensitivity}}, \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
Wynn = function(sensF, tol, maxIter=1e4) {
    tag = list(Wynn=list(sensF=sensF, tol=tol, maxIter=maxIter))

    defaults = attr(sensF, 'defaults')
    x = defaults[['x']]
    if (!isTRUE(all.equal(x, defaults[['desx']])))
        stop('sensitivity defaults for x and desx shall be equal')

    if (nrow(x) == 0)
        return(design(x, numeric(0), tag=tag))

    n = nrow(x)
    w = rep(1/n, n)
    tolBreak = F

    for (iIter in seq1(1, maxIter)) {
        sens = sensF(desw=w)

        maxIdx = which(sens == max(sens))
        if (length(maxIdx) != 1)
            maxIdx = sample(maxIdx, 1)

        dw = 1 / (iIter + 1)
        w = w * (1 - dw)
        w[maxIdx] = 0
        w[maxIdx] = 1 - sum(w) # equal to 'w[maxIdx] + dw'

        mSens = sens[maxIdx]
        if (mSens <= tol) {
            tolBreak = T
            break
        }
    }

    names(sens) = NULL
    tag$Wynn$tolBreak = tolBreak
    tag$Wynn$nIter = iIter
    return(design(x, w, tag=tag))
}


wPoint = function(x, w) {
    ## x = row matrix
    ## w = vector
    ## nrow(x) == length(w)
    return( apply(sweep(x, 1, w, '*'), 2, sum) / sum(w) )
}

#' Reduce Design
#'
#' \code{reduce} drops insignificant points and merges points in a certain neighbourhood.
#'
#' @param des a design.
#' @param distMax maximum euclidean distance between points to be merged.
#' @param wMin minimum weight a point shall have to be considered significant.
#'
#' @return \code{reduce} returns an object of \code{class} \code{"desigh"}.
#' See \code{\link{design}} for its structural definition.
#'
#' @seealso \code{\link{design}}
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

    tag = list(reduce=list(des=des, distMax=distMax, wMin=wMin))
    return(design(rx, rw, tag=tag))
}


#' Get Fisher Information
#'
#' \code{getM} returns the Fisher information corresponding to a model and a design.
#'
#' @param mod a model.
#' @param des a design.
#'
#' @return \code{getM} returns a named matrix, the Fisher information.
#'
#' @seealso \code{\link{param}}, \code{\link{design}}
#'
#' @examples ## see examples for param
#'
#' @export
getM = function(mod, des) {
    m = getm(mod, des$x)
    r = getM_(m, des$w)
    dimnames(r) = dimnames(mod$fisherI[[1]])
    return(r)
}


## from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
## Adds transparency to colours
add.alpha <- function(col, alpha=1){
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
}

getBaseDesign = function(des) {
    tag = des$tag
    if (!is.null(tag[['Wynn']]))
        return(des)
    if (!is.null(tag[['reduce']]))
        return(getBaseDesign(tag$reduce$des))
    return(NULL)
}

#' Plot Design
#'
#' \code{plot.desigh} creates a one-dimensional design plot, optionally together with a specified sensitivity curve.
#' If the design space has additional dimensions, the design is projected on a specified margin.
#'
#' @param x a design.
#' @param sensx (optional) a row matrix of points.
#' @param sens (optional) either a vector of sensitivities or a sensitivity function.
#' The latter shall rely on defaults, see \code{\link{Dsensitivity}} for details.
#' @param sensTol (optional) a single numeric.
#' Adds a horizontal line at this sensitivity level.
#' @param ... other arguments passed to plot.
#' @param margins a vector of indices, the dimensions to project on.
#' Defaults to \code{1}.
#' @param desSens if \code{TRUE} and \code{sens} is not specified then the sensitivity function which potentially was used in \code{Wynn} is taken as \code{sens}.
#' @param sensPch either a character vector of point 'characters' to add to the sensitivity curve or \code{NULL}.
#' @param sensArgs a list of arguments passed to draw calls related to the sensitivity.
#'
#' @references uses add.alpha from \url{http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html}
#'
#' @seealso \code{\link{design}}, \code{\link{Dsensitivity}}
#'
#' @examples ## see examples for param
#'
#' @export
plot.desigh = function(x, sensx=NULL, sens=NULL, sensTol=NULL, ..., margins=NULL, desSens=T, sensPch='+', sensArgs=list()) {
    points = function(..., axes) graphics::points(...)
    abline = function(..., axes) graphics::abline(...)
    axis = function(..., axes) graphics::axis(...)
    mtext = function(..., axes) graphics::mtext(...)

    des = x # workaround for S3 requirement

    args = list(...)

    if (is.null(margins))
        margins = 1
        #margins = 1:(ncol(des$model$x))
    if (1 < length(margins))
        stop('not yet implemented')

    ## marginal projection
    x = des$x
    w = des$w
    idcs = split(seq1(1, nrow(x)), lapply(margins, function(margin) x[, margin]), drop=T)
    x = x[sapply(idcs, function(i) i[1]), margins, drop=F]
    w = sapply(idcs, function(idcs) sum(w[idcs]))
    ord = roworder(x)
    x = x[ord,, drop=F]
    w = w[ord]

    ## lookup design sensF
    if (is.null(sens) && desSens) {
        d = getBaseDesign(des)
        sens = d$tag[['Wynn']]$sensF
        if (is.null(sensx) && !is.null(sens)) {
            sensx = attr(sens, 'defaults')[['x']]
            sens = sens(desw=d$w)
        }
    }

    ## prepare sensitivity
    if (!is.null(sens)) {
        if (is.function(sens)) {
            if (is.null(sensx)) {
                sensx = attr(sens, 'defaults')[['x']]
                sens = sens(desw=des$w, desx=des$x)
            } else
                sens = sens(x=sensx, desw=des$w, desx=des$x)
        }

        if (is.null(sensx))
            stop('if sens is a vector, sensx shall be specified')

        ## marginal projection
        idcs = split(seq1(1, nrow(sensx)), lapply(margins, function(margin) sensx[, margin]), drop=T)

        sensx = sensx[sapply(idcs, function(i) i[1]), margins, drop=F]
        sens = sapply(idcs, function(idcs) max(sens[idcs]))

        ord = roworder(sensx)
        sensx = sensx[ord,, drop=F]
        sens = sens[ord]

        ## add axis margin
        mar = par('mar')
        mar[4] = mar[2]
        par(mar=mar)
    } else {
        ## remove axis margin
        mar = par('mar')
        mar[4] = 2 + 0.1 # Warning, hard coded
        par(mar=mar)
    }

    if (length(margins) == 1) {
        xlab = colnames(x)
        if (is.null(xlab))
            xlab = paste('x[, c(', toString(margins), ')]', sep='')

        margs = modifyList(list(NA, xlim=range(x), ylim=c(0, 1), xlab=xlab, ylab='weight'), args)
        do.call(plot, margs)

        xlim = margs$xlim; ylim = margs$ylim

        if (!is.null(sens)) {
            par(new=T)

            dylim = range(sens)
            if (0 < dylim[1])
                dylim = c(0, dylim[2])
            else if (dylim[2] < 0)
                dylim = c(dylim[1], 0)
            if (isTRUE(sensTol < dylim[1]))
                dylim = c(sensTol, dylim[2])
            else if (isTRUE(dylim[2] < sensTol))
                dylim = c(dylim[1], sensTol)

            defaultYlab = ifelse(ncol(des$x) == 1, 'sensitivity', 'maximum sensitivity')
            margs = modifyList(list(sensx, sens, type='l', ylim=dylim, ylab=defaultYlab), sensArgs)
            ylab = margs$ylab
            margs = modifyList(margs, list(xlim=xlim, axes=F, xlab='', ylab=''))
            do.call(plot, margs)

            if (!is.null(sensTol)) {
                margs = modifyList(list(h=sensTol, col='black', lty=2), sensArgs)
                #margs$col = add.alpha(margs$col, 0.33)
                do.call(abline, margs)
            }

            if (!isTRUE(sensArgs[['axes']] == F)) {
                margs = modifyList(list(4), sensArgs)
                if (!is.null(margs[['col']]) && is.null(margs[['col.axis']]))
                    margs$col.axis = margs$col
                do.call(axis, margs)
            }

            margs = modifyList(list(ylab, side=4, line=3), sensArgs)
            do.call(mtext, margs)

            if (!(is.null(sensPch) || identical(sensPch, ''))) {
                alpha = des$w / max(des$w)
                idcs = which(1/256 < alpha)
                px = des$x[idcs,]
                p = approx(sensx, sens, px)

                margs = modifyList(list(p$x, p$y, pch=sensPch, col='black'), args)
                margs$col = add.alpha(margs$col, alpha[idcs])
                do.call(points, margs)
            }
        }

        par(new=T)
        margs = modifyList(list(x, w, xlim=xlim, ylim=ylim, type='h'), args)
        margs = modifyList(margs, list(axes=F, xlab='', ylab=''))
        do.call(plot, margs)
    }
}


#' D Efficiency
#'
#' \code{Defficiency} computes the D-, D_s or D_A-efficiency measure for a design with respect to a reference design.
#'
#' Indices supplied to argument \code{A} correspond to the subset of parameters defined by argument \code{parNames}.
#'
#' D efficiency is defined as
#' \deqn{\left(\frac{\left|M(\xi,\bar{\theta})\right|}{\left|M(\xi^{*},\bar{\theta})\right|}\right)^{1/n}}{( det(M(\xi, \theta)) / det(M(\xi*, \theta)) )**(1/n)}
#' and D_A efficiency as
#' \deqn{\left(\frac{\left|A^{T}M(\xi^{*},\bar{\boldsymbol{\theta}})^{-1}A\right|}{\left|A^{T}M(\xi,\bar{\boldsymbol{\theta}})^{-1}A\right|}\right)^{1/s}}{( det(t(A) \%*\% solve(M(\xi*, \theta)) \%*\% A) / det(t(A) \%*\% solve(M(\xi, \theta)) \%*\% A) )**(1/s)}
#'
#' @param des a design.
#' @param ref a design, the reference.
#' @param mod a model.
#' @param A for \itemize{
#' \item D-efficiency: \code{NULL}
#' \item D_s-efficiency: a vector of names or indices, the subset of parameters of interest.
#' \item D_A-efficiency: either \itemize{
#'   \item directly: a matrix without row names.
#'   \item indirectly: a matrix with row names corresponding to the parameters.
#'   }
#' }
#' @param parNames a vector of names or indices, the subset of parameters to use.
#' Defaults to the parameters for which the Fisher information is available.
#'
#' @return \code{Defficiency} returns a single numeric.
#'
#' @seealso \code{\link{design}}, \code{\link{param}}
#'
#' @examples ## see examples for param
#'
#' @export
Defficiency = function(des, ref, mod, A=NULL, parNames=NULL) {
    tt = getDAPar(mod, A, parNames)
    parNames = tt$parNames
    A = tt$A

    m = getm(mod, des$x, parNames)
    M = getM_(m, des$w)

    m = getm(mod, ref$x, parNames)
    Mref = getM_(m, ref$w)

    if (is.null(A))
        return( (det(M) / det(Mref))**(1/nrow(M)) )

    return( (det(t(A) %*% solve(Mref) %*% A) / det(t(A) %*% solve(M) %*% A) )**(1/ncol(A)) )
}

