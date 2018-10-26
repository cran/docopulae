

#' Sequence Generation
#'
#' \code{seq1} is similar to \code{seq}, however \code{by} is strictly \code{1} by default and \code{integer(0)} is returned if the range is empty.
#'
#' @param from,to,by see \code{\link[base]{seq}}.
#'
#' @return \code{seq1} returns either \code{integer(0)} if range is empty or what an appropriate call to \code{seq} returns otherwise.
#'
#' See examples below.
#'
#' @example R/examples/seq1.R
#'
#' @seealso \code{\link[base]{seq}}
#'
#' @export
seq1 = function(from, to, by=1) {
    if (to == from)
        return(from)
    if (sign(to - from) != sign(by))
        return(integer(0))
    return(seq(from, to, by))
}


Derivf = function(f, x, ...) {
    ## x = named unique vector
    xx = mapply(function(n, v) { # split x into list of named values
        names(v) = n
        return(v)
    }, names(x), x, SIMPLIFY=F)

    r = rep(list(NULL), length(xx))
    for (i in seq1(1, length(xx))) {
        xi = xx[[i]]
        cat('d/d ', xi, '\n', sep='')
        r[[i]] = Deriv::Deriv(f, xi, ...)
    }
    names(r) = x
    return(r)
}

Deriv2f = function(f, x, ...) {
    ## x = named unique vector
    d = Derivf(f, x, ...)

    xx = mapply(function(n, v) { # split x into list of named values
        names(v) = n
        return(v)
    }, names(x), x, SIMPLIFY=F)
    names(xx) = x

    r = rep(list(NULL), length(x))
    names(r) = x
    r = replicate(length(x), r, simplify=F)
    names(r) = x

    for (i in seq1(1, length(d))) {
        xi = xx[[i]]
        cat('d^2/d ', xi, ' d ', xi, '\n', sep='')

        r[[xi]][[xi]] = Deriv::Deriv(d[[i]], xi, ...)
    }

    combs = combn(as.character(x), 2)
    for (i in seq1(1, ncol(combs))) {
        a = combs[1, i]
        b = combs[2, i]
        cat('d^2/d ', a, ' d ', b, '\n', sep='')

        n = Deriv::Deriv(d[[a]], xx[[b]], ...)
        r[[a]][[b]] = n
        r[[b]][[a]] = n
    }

    return(r)
}


mirrorMatrix = function(x) {
    ## transforms upper/lower diagonal matrix to full matrix
    r = x + t(x)
    diag(r) = diag(x)
    return(r)
}


## nested list helper
is_flat = function(x) !any(sapply(x, inherits, 'list'))

flatten = function(x) {
    if (!inherits(x, 'list'))
        return(list(x))
    if (is_flat(x))
        return(x)
    return(do.call(c, lapply(x, flatten)))
}

rlapply = function(X, FUN, ...) {
    isList = vapply(X, inherits, T, 'list')
    r = rep(list(NULL), length(X))
    r[isList] = lapply(X[isList], rlapply, FUN, ...)
    r[!isList] = lapply(X[!isList], FUN, ...)
    return(r)
}


zmin = function(x) if (length(x) == 0) 0 else min(x)
zmax = function(x) if (length(x) == 0) 0 else max(x)


lproduct = function(x) {
    ## product expands a list of lists
    if (length(x) == 0)
        return(list())
    idcs = as.matrix(do.call(expand.grid, lapply(x, seq))) # row matrix of idx combinations
    m = suppressWarnings(do.call(cbind, x)) # matrix of objects
    names = suppressWarnings(do.call(cbind, lapply(x, base::names))) # matrix of names
    off = (0:(ncol(m) - 1)) * nrow(m) # column offsets
    if (is.null(names))
        r = apply(idcs, 1, function(idcs) m[idcs + off])
    else
        r = apply(idcs, 1, function(idcs) {
            idcs = idcs + off
            r = m[idcs]
            base::names(r) = names[idcs]
            return(r)
        })
    return(r)
}


#' Integrate Alternative
#'
#' \code{integrateA} is a tolerance wrapper for \code{integrate}.
#' It allows \code{integrate} to reach the maximum number of subdivisions.
#'
#' See \code{\link[stats]{integrate}}.
#'
#' @param f,lower,upper,...,subdivisions,rel.tol,abs.tol,stop.on.error,keep.xy,aux see \code{\link[stats]{integrate}}.
#'
#' @seealso \code{\link[stats]{integrate}}
#'
#' @example R/examples/integrateA.R
#'
#' @export
integrateA = function(f, lower, upper, ..., subdivisions=100L, rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, stop.on.error=TRUE, keep.xy=FALSE, aux=NULL) {
    r = stats::integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
    if ( !(r[['message']] %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) {
            stop(r$message)
        }
    }
    return(r)
}


clusterDist = function(distMat, maxDist) {
    ## distMat = square matrix of pairwise distances between objects
    ## maxDist = maximum distance within a cluster
    ##
    ## in ascending order of rows/columns
    ## iteratively assigns objects to the nearest cluster
    ##
    ## returns a vector of cluster ids

    r = integer(nrow(distMat))
    if (length(r) == 0)
        return(r)

    r[1] = 1
    rr = 2

    for (i in seq1(2, ncol(distMat))) {
        ds = distMat[1:(i - 1), i] # distances to other objects
        j = which.min(ds)
        if (ds[j] <= maxDist)
            r[i] = r[j]
        else {
            r[i] = rr
            rr = rr + 1
        }
    }

    return(r)
}

clusterPeak = function(x, y, maxDist) {
    ## x = row matrix of points
    ## y = corresponding vector of values (of length nrow(x))
    ##
    ## in descending order of magnitude of y
    ## iteratively assigns points to the nearest cluster
    ##
    ## returns a vector of cluster ids

    ord = order(y, decreasing=T)
    x = x[ord,, drop=F]
    dists = as.matrix(dist(x))
    r = clusterDist(dists, maxDist)
    r = r[order(ord)]
    return(r)
}


#' Matrix Ordering Permutation
#'
#' \code{roworder} returns a permutation which rearranges the rows of its first argument into ascending order.
#'
#' @param x a matrix.
#' @param ... other arguments passed to \code{order}.
#'
#' @return \code{roworder} returns an integer vector.
#'
#' @seealso \code{\link[base]{order}}
#'
#' @example R/examples/roworder.R
#'
#' @export
roworder = function(x, ...) {
    cols = lapply(seq1(1, ncol(x)), function(i) x[,i])
    args = c(cols, ...)
    r = do.call(order, args)
    if (is.null(r)) {
        r = integer(0)
    }
    return(r)
}


#' Row Matching
#'
#' \code{rowmatch} returns a vector of the positions of (first) matches of the rows of its first argument in the rows of its second.
#'
#' \code{rowmatch} uses compiled C-code.
#'
#' @param x a row matrix of doubles, the rows to be matched.
#' @param table a row matrix of doubles, the rows to be matched against.
#' @param nomatch the value to be returned in the case when no match is found.
#' Note that it is coerced to \code{integer}.
#'
#' @return \code{rowmatch} returns an integer vector giving the position of the matching row in \code{table} for each row in \code{x}. And \code{nomatch} if there is no matching row.
#'
#' @seealso \code{\link[base]{match}}
#'
#' @example R/examples/rowmatch.R
#'
#' @export
#' @useDynLib docopulae, rowmatch_double, .registration = TRUE 
rowmatch = function(x, table, nomatch=NA_integer_) {
    if (!(is.double(x) && is.matrix(x) && is.double(table) && is.matrix(table)))
        stop('x and table shall be matrices of doubles')

    if (ncol(x) != ncol(table))
        stop('x and table shall have the same number of columns')

    nomatch = as.integer(nomatch)

    ordx = roworder(x)
    ordt = roworder(table)
    sx = x[ordx,, drop=F]
    st = table[ordt,, drop=F]

    i = integer(nrow(x))
    i = .C('rowmatch_double', sx, as.integer(nrow(sx)), as.integer(ncol(sx)), st, as.integer(nrow(st)), i, PACKAGE='docopulae')[[6]] + 1

    r = ordt[i][order(ordx)]
    r[is.na(r)] = nomatch
    return(r)
}


#' Determine Duplicate Rows
#'
#' \code{rowsduplicated} determines which rows of a matrix are duplicates of rows with smaller subscripts, and returns a logical vector indicating which rows are duplicates.
#'
#' \code{rowsduplicated} uses compiled C-code.
#'
#' @param x a row matrix of doubles.
#'
#' @return \code{rowsduplicated} returns a logical vector with one element for each row.
#'
#' @seealso \code{\link[base]{duplicated}}
#'
#' @export
#' @useDynLib docopulae, rowsduplicated_double, .registration = TRUE
rowsduplicated = function(x) {
    if (!(is.double(x) && is.matrix(x)))
        stop('x shall be a matrix of doubles')

    order_ = roworder(x)
    x = x[order_,, drop=F]

    r = logical(nrow(x))
    r = .C('rowsduplicated_double', x, as.integer(nrow(x)), as.integer(ncol(x)), r, PACKAGE='docopulae')[[4]]

    r = r[order(order_)]
    return(r)
}


#' Grow Grid
#'
#' \code{grow.grid} creates a data frame like \code{expand.grid}.
#' The order of rows is adjusted to represent a growing grid with respect to resolution.
#'
#' @param x a list of vectors.
#' @param random \code{TRUE} if order of rows within each level of resolution should be random.
#'
#' @return \code{grow.grid} returns a data frame like \code{expand.grid}.
#'
#' @seealso \code{\link{update.param}}
#'
#' @export
grow.grid = function(x, random=T) {
    if (length(x) == 0) {
        return(data.frame())
    }
    stages = lapply(x, get.stages)
    # ensure same length
    n = max(sapply(stages, length))
    stages = lapply(stages, function(stages) c(stages, rep(stages[length(stages)], n - length(stages))))
    #
    grids = do.call(function(...) mapply(expand.grid, ..., SIMPLIFY=F), stages)
    grids = lapply(grids, as.data.frame) # ensure data frames
    if (random) {
        grids = lapply(grids, function(grid) grid[sample.int(nrow(grid)),, drop=F])
    }
    r = unique(do.call(rbind, grids))
    rownames(r) = NULL
    return(r)
}

get.stages = function(x) {
    if (length(x) == 0) {
        return(list())
    }
    n = ceiling(log2(length(x))*2 + 2) # overestimate number of stages
    r = rep(list(NULL), n)
    step = length(x) - 1 # start with lowest resolution
    for (i in 1:length(r)) {
        j = unique(round(seq(1, length(x), step))) # get indices
        r[[i]] = x[j]
        if (length(j) == length(x)) {
            break
        }
        step = step / 2 # increase resolution
    }
    r = r[1:i] # remove empty
    return(r)
}

