## ideas
## - implement simplify/merge


#' Dimension Type Attribute Values
#'
#' A dimension object is identified by its dimension type attribute \code{"nint_dtype"}.
#' On creation it is set to one of the following.
#' See dimension types in "See Also" below.
#'
#' @format integer
#' @seealso \code{\link{nint_scatDim}}, \code{\link{nint_gridDim}}, \code{\link{nint_intvDim}}, \code{\link{nint_funcDim}}, \code{\link{nint_space}}
#'
#' @name nint_TYPE
NULL

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_SCAT_DIM # = 1
#'
#' @export
nint_TYPE_SCAT_DIM = 1

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_GRID_DIM # = 2
#'
#' @export
nint_TYPE_GRID_DIM = 2

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_INTV_DIM # = 3
#'
#' @export
nint_TYPE_INTV_DIM = 3

#' @rdname nint_TYPE
#'
#' @usage nint_TYPE_FUNC_DIM # = 4
#'
#' @export
nint_TYPE_FUNC_DIM = 4


#' Scatter Dimension
#'
#' \code{nint_scatDim} is defined by a sequence of values.
#' Together with other scatter dimensions it defines a sparse grid.
#'
#' Imagine using \code{cbind} to create a row matrix of points.
#'
#' @param x a vector of any type.
#'
#' @return \code{nint_scatDim} returns its argument with the dimension type attribute set to \code{nint_TYPE_SCAT_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_scatDim = function(x) {
    x = as.vector(x)
    attr(x, 'nint_dtype') = nint_TYPE_SCAT_DIM
    return(x)
}

#' Grid Dimension
#'
#' \code{nint_gridDim} is defined by a sequence of values.
#' Together with other grid dimensions it defines a dense grid.
#'
#' Imagine using \code{expand.grid} to create a row matrix of points.
#'
#' @param x a vector of any type.
#'
#' @return \code{nint_scatDim} returns its argument with the dimension type attribute set to \code{nint_TYPE_GRID_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_gridDim = function(x) {
    x = as.vector(x)
    attr(x, 'nint_dtype') = nint_TYPE_GRID_DIM
    return(x)
}

#' Interval Dimension
#'
#' \code{nint_intvDim} defines a fixed interval.
#' The bounds may be (negative) \code{Inf}.
#'
#' @param x either a single numeric, the lower bound, or a vector of length 2, the lower and upper bound.
#' @param b the upper bound if \code{x} is the lower bound.
#'
#' @return \code{nint_intvDim} returns a vector of length 2 with the dimension type attribute set to \code{nint_TYPE_INTV_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_intvDim = function(x, b=NULL) {
    if (is.null(b))
        x = as.vector(x)
    else
        x = c(x, b)
    if (length(x) != 2) stop('exactly two values shall be specified as bounds')
    attr(x, 'nint_dtype') = nint_TYPE_INTV_DIM
    return(x)
}

#' Function Dimension
#'
#' \code{nint_funcDim} defines a functionally dependent dimension.
#' It shall depend solely on the previous dimensions.
#'
#' Obviously if \code{x} returns an object of type \code{nint_intvDim} the dimension is continuous, and discrete otherwise.
#'
#' As the argument to \code{x} is only partially defined the user has to ensure that the function solely depends on values up to the current dimension.
#'
#' @param x \code{function(x)}, where \code{x} is the partially realized point in the space.
#' It shall return an object of type \code{nint_intvDim} or a vector.
#'
#' @return \code{nint_scatDim} returns its argument with the dimension type attribute set to \code{nint_TYPE_FUNC_DIM}.
#'
#' @seealso \code{\link{nint_TYPE}}, \code{\link{nint_space}}
#'
#' @export
nint_funcDim = function(x) {
    x = as.function(x)
    attr(x, 'nint_dtype') = nint_TYPE_FUNC_DIM
    return(x)
}


nint_dtype = function(x) attr(x, 'nint_dtype')

nint_toString_scatDim = function(x) paste('s(', toString(x), ')', sep='')
nint_toString_gridDim = toString
nint_toString_intvDim = function(x) paste('[', x[1], ', ', x[2], ']', sep='')
nint_toString_funcDim = function(x) paste(capture.output(print(attr(x, 'srcref'))), collapse='\n')

nint_toString_dim = function(x) {
    type = nint_dtype(x)
    if (is.null(type))
        return('<unknown dim type>')
    return(switch(type, nint_toString_scatDim, nint_toString_gridDim, nint_toString_intvDim, nint_toString_funcDim)(x))
}

nint_print_spaceDims = function(x) {
    x = flatten(x)
    x = lapply(x, nint_toString_dim)
    x = paste(x, collapse=' or ')
    cat(x, '\n')
}


#' Space Validation Errors
#'
#' Error codes for space validation.
#'
#' @format integer
#' @seealso \code{\link{nint_validateSpace}}, \code{\link{nint_space}}
#'
#' @name nint_ERROR
NULL

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_DIM_TYPE # = -1001
#'
#' @details \code{nint_ERROR_DIM_TYPE}: dimension type attribute does not exist or is invalid.
nint_ERROR_DIM_TYPE = -1001

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_SCATTER_LENGTH # = -1002
#'
#' @details \code{nint_ERROR_SCATTER_LENGTH}: scatter dimensions have different lengths.
nint_ERROR_SCATTER_LENGTH = -1002

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_SPACE_TYPE # = -1003
#'
#' @details \code{nint_ERROR_SPACE_TYPE}: object not of type \code{"nint_space"}.
nint_ERROR_SPACE_TYPE = -1003

#' @rdname nint_ERROR
#'
#' @usage nint_ERROR_SPACE_DIM # = -1004
#'
#' @details \code{nint_ERROR_SPACE_DIM}: subspaces have different number of dimensions.
nint_ERROR_SPACE_DIM = -1004


nint_validateSpaceDim = function(x, refl) {
    ## check - dim type
    ##       - scat length
    type = nint_dtype(x)
    if (is.null(type))
        return(nint_ERROR_DIM_TYPE)
    pass = function(x, refl) 0
    f = switch(type, function(x, refl) ifelse(length(x) != refl, nint_ERROR_SCATTER_LENGTH, 0), pass, pass, pass)
    if (is.null(f))
        return(nint_ERROR_DIM_TYPE)
    return(f(x, refl))
}

nint_validateSpaceDims = function(x, refl) zmin(sapply(flatten(x), nint_validateSpaceDim, refl))


#' Space
#'
#' \code{nint_space} defines an n-dimensional space as a list of dimensions.
#' A space may consist of subspaces.
#' A space without subspaces is called true subspace.
#'
#' If a space contains at least one list structure of dimension objects it consists of subspaces.
#' Each subspace is then defined by a combination of dimension objects along the dimensions.
#' See \code{\link{nint_expandSpace}} on how to expand a space to true subspaces.
#'
#' @param ... dimensions each of which may be an actual dimension object or a list structure of dimension objects.
#'
#' @return \code{nint_space} returns an object of \code{class} \code{"nint_space"}.
#' An object of \code{class} \code{"nint_space"} is an ordered list of dimension objects.
#'
#' @seealso \code{\link{nint_scatDim}}, \code{\link{nint_gridDim}}, \code{\link{nint_intvDim}}, \code{\link{nint_funcDim}}, \code{\link{nint_integrate}}, \code{\link{nint_validateSpace}}, \code{\link{nint_expandSpace}}, \code{\link{fisherI}}
#'
#' @example R/examples/nint_space.R
#'
#' @export
nint_space = function(...) {
    r = list(...)
    class(r) = 'nint_space'
    return(r)
}

#' Print Space
#'
#' \code{print.nint_space} prints a space in a convenient way.
#'
#' Each line represents a dimension.
#' Format: "<dim idx>: <dim repr>".
#' Each dimension has its own representation which should be easy to understand.
#' \code{nint_scatDim} representations are marked by \code{"s()"}.
#'
#' @param x a space.
#' @param ... ignored.
#'
#' @seealso \code{\link{nint_space}}
#'
#' @export
print.nint_space = function(x, ...) {
    for (i in seq1(1, length(x))) {
        cat(i, ': ', sep='')
        nint_print_spaceDims(x[[i]])
    }
}

#' Validate Space
#'
#' \code{nint_validateSpace} performs a couple of checks on a space or list structure of spaces to ensure it is properly defined.
#'
#' @param x a space or list structure of spaces.
#'
#' @return \code{nint_validateSpace} returns 0 if everything is fine, or an error code.
#' See \code{\link{nint_ERROR}}.
#'
#' @seealso \code{\link{nint_ERROR}}, \code{\link{nint_space}}
#'
#' @example R/examples/nint_validateSpace.R
#'
#' @export
nint_validateSpace = function(x) {
    x = flatten(x)
    if (length(x) == 0)
        return(0)
    return(min(sapply(x, nint_validateSpace_, x[[1]])))
}

nint_validateSpace_ = function(x, refs=NULL) {
    if (!inherits(x, 'nint_space')) # check class
        return(nint_ERROR_SPACE_TYPE)
    if (length(x) != length(refs)) # check nint_space dimension
        return(nint_ERROR_SPACE_DIM)
    x = lapply(x, flatten)
    return(zmin( sapply(x, nint_validateSpaceDims, nint_validateSpace_getRefl(x)) ))
}

nint_validateSpace_getRefl = function(x) {
    return(zmax(sapply(x, function(x) zmax(sapply(x, nint_validateSpace_getRefl_) )) ))
}

nint_validateSpace_getRefl_ = function(x) {
    type = nint_dtype(x)
    if (is.null(type) || type != nint_TYPE_SCAT_DIM)
        return(0)
    return(length(x))
}

#' Expand Space
#'
#' \code{nint_expandSpace} expands a space or list structure of spaces to a list of true subspaces.
#'
#' @param x a space or list structure of spaces.
#'
#' @return \code{nint_expandSpace} returns a list of spaces.
#' Each space is a true subspace.
#'
#' @seealso \code{\link{nint_space}}
#'
#' @example R/examples/nint_expandSpace.R
#'
#' @export
nint_expandSpace = function(x) {
    return(flatten(lapply(flatten(x), nint_expandSpace_) ))
}

nint_expandSpace_ = function(x) {
    return( lapply(lproduct(lapply(x, flatten)), function(x) do.call(nint_space, x)) )
}


## Create Integration Space
##
## \code{nint_ispace} transforms a true subspace into a data structure which is used by \code{nint_integrate} to efficiently integrate over the entire space.
##
## Assumes interchangeability of dimensions (except function dimensions).
##
## @param x a true subspace.
##
## @return \code{nint_ispace} returns a list.
##
## @seealso \code{\link{nint_space}}, \code{\link{nint_integrate}}
##
## @example R/examples/nint_ispace.R
##
## @export
nint_ispace = function(x) {
    ## - group dimensions by type and put them in the following order
    ##   - scatter (s)
    ##   - grid (g)
    ##   - interval (i)
    ##   - function (f)
    ## - bind scattered dimensions by column
    ## - expand grid dimensions to grid
    ## - bind interval bounds by row
    ## - create function list
    ##
    ## - result := list(type=list(i=idcs, g=data))
    ##   type := c('s', 'g', 'i', 'f')
    ##   data := matrix      if type == 's'
    ##        or data.frame  if type == 'g' or type == 'i'
    ##        or list        if type == 'f'

    if (length(x) == 0)
        return(list())

    rr = list()
    si = integer(0) # scatter
    gi = integer(0) # grid
    ii = integer(0) # interval
    fi = integer(0) # function

    for (i in seq1(1, length(x))) {
        xx = x[[i]]
        if (inherits(xx, 'list'))
            stop('argument is no true subspace')
        type = nint_dtype(xx)
        if (is.null(type))
            stop('unknown dimension type')

        if (type == nint_TYPE_FUNC_DIM) {
            if (length(si) != 0) rr = c(rr, list(s=si))
            if (length(gi) != 0) rr = c(rr, list(g=gi))
            if (length(ii) != 0) rr = c(rr, list(i=ii))
            si = integer(0)
            gi = integer(0)
            ii = integer(0)
            fi = c(fi, i)
            next
        }
        if (length(fi) != 0) {
            rr = c(rr, list(f=fi))
            fi = integer(0)
        }

        if (type == nint_TYPE_SCAT_DIM) {
            si = c(si, i)
        } else if (type == nint_TYPE_GRID_DIM) {
            gi = c(gi, i)
        } else if (type == nint_TYPE_INTV_DIM) {
            ii = c(ii, i)
        } else {
            stop('unknown dimension type')
        }
    }

    if (length(si) != 0) rr = c(rr, list(s=si))
    if (length(gi) != 0) rr = c(rr, list(g=gi))
    if (length(ii) != 0) rr = c(rr, list(i=ii))
    if (length(fi) != 0) rr = c(rr, list(f=fi))

    ts = function(i) do.call(cbind, x[i])
    tg = function(i) do.call(expand.grid, x[i])
    ti = function(i) do.call(rbind, x[i])
    tf = function(i) x[i]
    tt = list(s=ts, g=tg, i=ti, f=tf)
    r = mapply(function(type, i) list(i=i, g=tt[[type]](i)), names(rr), rr, SIMPLIFY=F)

    ## create matrix containing the depth (first column) and position (second column) for each dimension
    #depth = rep.int(1:length(rr), sapply(rr, length))[order(unlist(rr))]
    #attr(r, 'depth') = cbind(depth, mapply(match, 1:length(x), lapply(depth, function(d) r[[d]]$i)), deparse.level=0)
    return(r)
}

nint_ispaces = function(x) {
    r = nint_validateSpace(x)
    if (r != 0)
        stop(r)
    x = nint_expandSpace(x)
    return(lapply(x, nint_ispace))
}


## not needed
## @export
#nint_applyToSpace = function(space, d, f) {
    #r = lapply(flatten(space), nint_applyToSpace_, d, f)
    #r = flatten(r)
    #if (length(r) == 1)
        #return(r[[1]])
    #return(r)
#}

#insert = function(x, idcs, l) {
    ## x .. list of objects
    ## idcs .. vector of indices
    ## l .. list to insert to
    #r = l
    #d = length(x) - length(idcs)
    #if (d == 0) {
        #r[idcs] = x
    #} else if (d < 0) { # more idcs than items
        #r[idcs[seq1(1, length(x))]] = x # set
        #r = r[-idcs[seq1(length(x) + 1, length(idcs))]] # del
        #r = do.call(nint_space, r)
    #} else { # more items than idcs
        #r[idcs] = x[seq1(1, length(idcs))] # set
        #last = idcs[length(idcs)]
        #r = c(r[seq1(1, last)], x[seq1(length(idcs) + 1, length(x))], r[seq1(last + 1, length(r))]) # insert
        #r = do.call(nint_space, r)
    #}
    #return(r)
#}

#nint_applyToSpace_ = function(space, d, f) {
    #dds = lapply(space[d], flatten) # list(list(dim))
    #r = lapply(lproduct(dds), f) # list(obj)
    #r = lapply(r, function(x) if (inherits(x, 'list')) x else list(x)) # list(list(dim))
    #r = lapply(r, insert, d, space) # list(space)
    #return(r)
#}


#' Tangent Transform
#'
#' \code{nint_tanTransform} creates the transformation \code{g(x) = atan((x - center)/scale)} to be used in \code{nint_transform}.
#'
#' @param center,scale see \code{g(x)}.
#' @param dIdcs an integer vector of indices, the dimensions to transform.
#'
#' @return \code{nint_tanTransform} returns a named list of two functions \code{"g"} and \code{"giDgi"} as required by \code{nint_transform}.
#'
#' @seealso \code{\link{nint_transform}}
#'
#' @example R/examples/nint_tanTransform.R
#'
#' @export
nint_tanTransform = function(center, scale, dIdcs=NULL) {
    tt = list(center, scale, dIdcs)
    g = function(x) atan((x - center)/scale)
    giDgi = function(y) {
        ty = tan(y)
        list(ty*scale + center,
             scale*(1 + ty^2)) # better than scale/(cos(y)^2)
    }
    r = list(g=g, giDgi=giDgi)
    if (!is.null(dIdcs))
        r$dIdcs = dIdcs
    return(r)
}


transform_error1 = 'transformed function dimensions shall return interval dimensions'

transform_funcDim = function(d, deps=NULL, g=NULL, xx=NULL, i=NULL, o=NULL, addNoDeps=F) {
    if (!is.null(g)) { # g
        tt = list(d, deps, g, xx, i, o)

        r = function(x) { # transform
            nx = x
            for (dep in deps)
                nx[ dep[['dIdcs']] ] = dep[['giD']](x[ dep[['inDIdcs']] ])[[1]]
            nd = d(nx)

            dtype = nint_dtype(nd)
            if (is.null(dtype) || dtype != nint_TYPE_INTV_DIM)
                stop(transform_error1)

            xx[i] = nd[1]
            a = g(xx)[o]
            xx[i] = nd[2]
            b = g(xx)[o]
            if (b < a)
                return(nint_intvDim(b, a))
            return(nint_intvDim(a, b))
        }

    } else if (!is.null(deps) && length(deps) != 0) { # !g, deps
        tt = list(d, deps)

        r = function(x) { # passive transform only
            nx = x
            for (dep in deps)
                nx[ dep[['dIdcs']] ] = dep[['giD']](x[ dep[['inDIdcs']] ])[[1]]
            return(d(nx))
        }

    } else # !g, !deps
        r = d

    r = nint_funcDim(r)
    if (addNoDeps) {
        if (!is.null(deps) && length(deps) != 0)
            attr(r, 'noDeps') = transform_funcDim(d, deps=NULL, g=g, xx=xx, i=i, o=o, addNoDeps=F)
        else
            attr(r, 'noDeps') = r
    }
    return(r)
}

transform_tranIOd = function(trans, maxDIdx) {
    dIdcs = lapply(trans, function(tran) tran$dIdcs)

    r = rep(list(NULL), maxDIdx)
    r[unlist(dIdcs, recursive=F)] = unlist(lapply(trans, function(tran)
        mapply(function(i, o)
            list(tran=tran, i=i, o=o)
        , match(tran$dIdcs, tran$inDIdcs), seq1(1, length(tran$dIdcs)), SIMPLIFY=F)
    ), recursive=F)

    return(r)
}

transform_space = function(space, trans, funcDimToF, defaultToF) {
    error1 = 'dimensions shall be either of type interval or function'

    dIdcs = lapply(trans, function(tran) tran$dIdcs)
    minDIdcs = sapply(dIdcs, min)
    maxDIdx = max(zmax(unlist(dIdcs, recursive=F)), funcDimToF)

    # per dimension
    depsd = lapply(seq1(1, maxDIdx), function(dIdx) trans[minDIdcs < dIdx]) # passive funcDim transformations
    tranIOd = transform_tranIOd(trans, maxDIdx)
    if (defaultToF)
        isToFd = rep(T, maxDIdx)
    else
        isToFd = seq1(1, maxDIdx) %in% funcDimToF

    r = rlapply(list(space), function(space) { # for each true space
        nDims = lapply(seq1(1, length(space)), function(dIdx) { # for each dimension
            d = space[[dIdx]]
            if (dIdx <= maxDIdx) {
                deps = depsd[[dIdx]]
                tranIO = tranIOd[[dIdx]]
                isToF = isToFd[dIdx]
            } else {
                deps = trans
                tranIO = NULL
                isToF = defaultToF
            }
            if (length(deps) == 0 && is.null(tranIO))
                return(d)

            if (!is.null(tranIO)) { # transform
                tran = tranIO$tran
                i = tranIO$i
                o = tranIO$o
                xx = rep(NA_real_, length(tran$inDIdcs))

                r = rlapply(list(d), function(d) { # for each true dimension
                    dtype = nint_dtype(d)
                    if (is.null(dtype))
                        stop(error1)

                    if (dtype == nint_TYPE_INTV_DIM) { # transform interval dimension
                        xx[i] = d[1]
                        a = tran$g(xx)[o]
                        xx[i] = d[2]
                        b = tran$g(xx)[o]
                        if (b < a)
                            return(nint_intvDim(b, a))
                        return(nint_intvDim(a, b))
                    }

                    if (dtype != nint_TYPE_FUNC_DIM)
                        stop(error1)
                    return(transform_funcDim(d, deps=deps, g=tran$g, xx=xx, i=i, o=o, addNoDeps=isToF)) # transform function dimension
                })[[1]]

            } else { # passive transform
                r = rlapply(list(d), function(d) { # for each true dimension
                    dtype = nint_dtype(d)
                    if (is.null(dtype) || dtype != nint_TYPE_FUNC_DIM)
                        return(d)

                    return(transform_funcDim(d, deps=deps, addNoDeps=isToF))
                })[[1]]
            }
            return(r)
        })

        r = space
        r[] = nDims
        return(r)
    })[[1]]

    return(r)
}

transform_f = function(f, trans, zeroInf) {
    tt = list(f, trans, zeroInf)

    trans = lapply(trans, function(tran) {
        tran$valid = which(!is.na(tran$dIdcs))
        return(tran)
    })

    r = function(x, ...) {
        J = 1
        iJ = 1
        for (tran in trans) {
            tt = tran[['giD']](x[ tran[['inDIdcs']] ])
            valid = tran[['valid']]
            x[tran[['dIdcs']][valid]] = tt[[1]][valid]
            if (tran[['isDgi']])
                J = J * prod(tt[[2]], na.rm=T)
            else
                iJ = iJ * prod(tt[[2]], na.rm=T)
        }
        J = J / iJ

        v = f(x, ...)
        if (is.infinite(J))
            if (v == 0)
                return(zeroInf)
        return(v*abs(J))
    }
    return(r)
}

transform_funcDimToF = function(d, dIdx, deps=NULL) {
    tt = list(d, dIdx, deps)

    giDgi = function(y) {
        ny = y
        for (dep in deps)
            ny[ dep[['dIdcs']] ] = dep[['giD']](y[ dep[['inDIdcs']] ])[[1]]
        nd = d(ny)

        dtype = nint_dtype(nd)
        if (is.null(dtype) || dtype != nint_TYPE_INTV_DIM)
            stop(transform_error1)

        t1 = nd[2] - nd[1]
        return(list(y[dIdx]*t1 + nd[1],
                    t1))
    }
    return(giDgi)
}

#' Transform Integral
#'
#' \code{nint_transform} applies monotonic transformations to an integrand and a space or list structure of spaces.
#' Common use cases include the probability integral transform, the transformation of infinite limits to finite ones and function dimensions to interval dimensions.
#'
#' Interval dimensions and function dimensions returning interval dimensions only.
#'
#' If a transformation is vector valued, that is \code{y = c(y1, ..., yn) = g(c(x1, ..., xn))}, then each component of \code{y} shall exclusively depend on the corresponding component of \code{x}.
#' So \code{y[i] = g[i](x[i])} for an implicit function \code{g[i]}.
#'
#' The transformation of function dimensions to interval dimensions is performed after the transformations defined by \code{trans}.
#' Consecutive linear transformations, \code{g(x[dIdx]) = (x[dIdx] - d(x)[1])/(d(x)[2] - d(x)[1])} where \code{d} is the function dimension at dimension \code{dIdx}, are used.
#' Deciding against this transformation probably leads to considerable loss in computational performance.
#'
#' @param f \code{function(x, ...)}, an integrand.
#' @param space a space or list structure of spaces.
#' @param trans a list of named lists, each containing \code{dIdcs}, \code{g} and \code{giDgi} or \code{giDg}, where\itemize{
#' \item \code{dIdcs} is an integer vector of indices, the dimensions to transform
#' \item \code{g=function(x[dIdcs])} mapping \code{x[dIdcs]} to \code{y}
#' \item \code{giDgi=function(y)} returning a list of two, the inverse \code{gi(y) = x[dIdcs]} and the first derivatives of \code{gi(y)} with respect to \code{y}
#' \item or \code{giDg=function(y)} returning the inverse and the first derivatives of \code{g(x[dIdcs])} with respect to \code{x[dIdcs]}.
#' }
#' @param funcDimToF an integer vector of indices, the dimensions to look for function dimensions to transform to interval dimensions.
#' \code{0} indicates all dimensions.
#' @param zeroInf a single value, used when \code{f} returns \code{0} and the Jacobian is infinite.
#'
#' @return \code{nint_transform} returns either a named list containing the transformed integrand and space, or a list of such.
#'
#' @seealso \code{\link{nint_integrate}}, \code{\link{nint_space}}, \code{\link{nint_tanTransform}}, \code{\link{fisherI}}
#'
#' @example R/examples/nint_transform.R
#'
#' @export
nint_transform = function(f, space, trans, funcDimToF=0, zeroInf=0) {
    tt = list(f, space, trans, funcDimToF, zeroInf)

    # prepare
    trans = lapply(trans, function(tran) {
        if (!is.null(tran[['giDgi']])) {
            isDgi = T
            giD = tran$giDgi
        } else if (!is.null(tran[['giDg']])) {
            isDgi = F
            giD = tran$giDg
        } else
            stop('transformations shall contain giDg or giDgi')
        dIdcs = tran[['dIdcs']]
        list(inDIdcs=dIdcs, g=tran[['g']], giD=giD, dIdcs=dIdcs, isDgi=isDgi)
    })

    dIdcs = lapply(trans, function(tran) tran$dIdcs)
    if (any(duplicated(unlist(dIdcs, recursive=F))))
        stop('each dimension shall be specified only once')

    defaultToF = isTRUE(funcDimToF == 0)

    if (length(funcDimToF) == 0) {
        # transform f
        rf = transform_f(f, trans, zeroInf)

        # transform space
        rspace = space

        if (length(trans) != 0)
            rspace = transform_space(rspace, trans, funcDimToF, defaultToF)

    } else {
        funcDimToF = sort(funcDimToF)
        nDim = nint_intvDim(0, 1)
        g = function(x) stop('shouldn\'t happen')

        maxDIdx = max(zmax(unlist(dIdcs, recursive=F)), funcDimToF)
        tranIOd = transform_tranIOd(trans, maxDIdx)

        # prepare space
        rspace = space

        if (length(trans) != 0)
            rspace = transform_space(rspace, trans, funcDimToF, defaultToF)

        r = lapply(flatten(rspace), function(space) { # for each true space
            if (defaultToF)
                funcDimToF = seq1(1, length(space))

            # split funcDim from others
            dims = lapply(space[funcDimToF], function(d) { # for each dimension
                dd = flatten(d)

                isFunc2 = sapply(dd, function(d) {
                    dtype = nint_dtype(d)
                    return(!is.null(dtype) && dtype == nint_TYPE_FUNC_DIM)
                })
                if (all(isFunc2)) {
                    names(dd) = rep('f', length(dd))
                    return(dd)
                }

                others = dd[!isFunc2]
                if (length(others) == 1)
                    others = others[[1]]

                dd = dd[isFunc2]
                names(dd) = rep('f', length(dd))

                return(c(dd, list(o=others)))
            })

            lapply(lproduct(dims), function(dims) { # for each combination of dimensions
                isFunc2 = (names(dims) == 'f')
                funcDimToF2 = funcDimToF[isFunc2]
                dims2 = dims[isFunc2]

                # create trans
                preTrans = lapply(trans, function(tran) {
                    i = which(tran$dIdcs %in% funcDimToF2)
                    if (length(i) == length(tran$dIdcs))
                        return(NULL)
                    tran$inDIdcs[i] = NA
                    tran$dIdcs[i] = NA
                    return(tran)
                })
                preTrans = preTrans[!sapply(preTrans, is.null)]

                trans1 = list()
                for (i in seq1(1, length(funcDimToF2))) {
                    dIdx = funcDimToF2[i]
                    giDgi = transform_funcDimToF(dims2[[i]], dIdx, deps=trans1)
                    n = list(inDIdcs=1:dIdx, g=g, giD=giDgi, dIdcs=dIdx, isDgi=T)
                    trans1 = c(trans1, list(n))
                }

                trans2 = mapply(function(d, dIdx) {
                    noDeps = attr(d, 'noDeps')
                    if (is.null(noDeps))
                        noDeps = d
                    giDgi = transform_funcDimToF(noDeps, dIdx)
                    list(inDIdcs=1:dIdx, g=g, giD=giDgi, dIdcs=dIdx, isDgi=T)
                }, dims2, funcDimToF2, SIMPLIFY=F)

                postTrans = lapply(funcDimToF2, function(dIdx) {
                    if (maxDIdx < dIdx)
                        return(NULL)
                    tranIO = tranIOd[[dIdx]]
                    if (is.null(tranIO))
                        return(NULL)
                    tran = tranIO$tran
                    tran$inDIdcs[-tranIO$i] = NA
                    tran$dIdcs[-tranIO$o] = NA
                    return(tran)
                })

                interleaved = rep(list(NULL), length(trans2) + length(postTrans))
                interleaved[seq(1, by=2, length.out=length(trans2))] = trans2
                interleaved[seq(2, by=2, length.out=length(postTrans))] = postTrans
                interleaved = interleaved[!sapply(interleaved, is.null)]

                # transform f
                rf2 = transform_f(f, c(preTrans, interleaved), zeroInf)

                # transform space
                rspace2 = space
                rspace2[funcDimToF] = dims
                rspace2 = transform_space(rspace2, trans1, c(), F)
                rspace2[funcDimToF2] = rep(list(nDim), length(funcDimToF2))

                list(f=rf2, space=rspace2)
            })
        })

        r = unlist(r, recursive=F)
        if (length(r) == 1)
            r = r[[1]]
        return(r)
    }

    return(list(f=rf, space=rspace))
}


#' Integrate Hypercube
#'
#' Interface to the integration over interval dimensions.
#'
#' \code{nint_integrate} uses \code{nint_integrateNCube} to handle interval dimensions.
#' See examples below on how to deploy different solutions.
#'
#' @usage nint_integrateNCube(f, lowerLimit, upperLimit, ...)
#'
#' @param f the scalar-valued wrapper function to be integrated.
#' @param lowerLimit the lower limits of integration.
#' @param upperLimit the upper limits of integration.
#' @param ... other arguments passed to \code{f}.
#'
#' @return \code{nint_integrateNCube} returns a single numeric.
#'
#' @seealso \code{\link{nint_integrate}}
#'
#' @example R/examples/nint_integrateNCube.R
#'
#' @name nint_integrateNCube
#'
#' @export
NULL

#' @details The function built by \code{nint_integrateNCube_integrate} calls \code{integrate} (argument) recursively.
#' The number of function evaluations therefore increases exponentially with the number of dimensions (\code{(subdivisions * 21) ** D} if \code{integrate}, the default, is used).
#' At the moment it is the default method because no additional package is required.
#' However, you most likely want to consider different solutions.
#'
#' @param integrate \code{function(f, lowerLimit, upperLimit, ...)} which calls \code{integrate}.
#'
#' @return \code{nint_integrateNCube_integrate} returns a recursive implementation for \code{nint_integrateNCube} based on one dimensional integration.
#'
#' @seealso \code{\link{integrateA}}, \code{\link[stats]{integrate}}
#'
#' @rdname nint_integrateNCube
#'
#' @export
nint_integrateNCube_integrate = function(integrate) {
    tt = integrate
    r = function(f, lowerLimit, upperLimit, ...) {
        maxDepth = length(lowerLimit)

        x = rep(0, maxDepth)
        d = 0
        g = Vectorize(function(xx, ...) {
            x[d] <<- xx
            if (d == maxDepth) {
                return(f(x, ...))
            }

            d <<- d + 1
            r = integrate(g, lowerLimit[d], upperLimit[d], ...)[['value']]
            d <<- d - 1
            return(r)
        }, 'xx')

        return(g(0, ...))
    }
    return(r)
}

#' @details The function built by \code{nint_integrateNCube_cubature} is a trivial wrapper for \code{cubature::adaptIntegrate}.
#'
#' @param adaptIntegrate \code{function(f, lowerLimit, upperLimit, ...)} which calls \code{cubature::adaptIntegrate}.
#'
#' @return \code{nint_integrateNCube_cubature} returns a trivial implementation for \code{nint_integrateNCube} indirectly based on \code{cubature::adaptIntegrate}.
#'
#' @seealso \code{\link[cubature]{adaptIntegrate}} in package \pkg{cubature}
#'
#' @rdname nint_integrateNCube
#'
#' @export
nint_integrateNCube_cubature = function(adaptIntegrate) {
    tt = adaptIntegrate
    r = function(f, lowerLimit, upperLimit, ...) {
        return(adaptIntegrate(f, lowerLimit, upperLimit, ...)[['integral']])
    }
    return(r)
}

#' @details The function built by \code{nint_integrateNCube_SparseGrid} is an almost trivial wrapper for \code{SparseGrid::createIntegrationGrid}.
#' It scales the grid to the integration region.
#'
#' @param createIntegrationGrid \code{function(dimension)} which calls \code{SparseGrid::createIntegrationGrid}.
#'
#' @return \code{nint_integrateNCube_SparseGrid} returns an implementation for \code{nint_integrateNCube} indirectly based on \code{SparseGrid::createIntegrationGrid}.
#'
#' @seealso \code{\link[SparseGrid]{createIntegrationGrid}} in package \pkg{SparseGrid}
#'
#' @rdname nint_integrateNCube
#'
#' @export
nint_integrateNCube_SparseGrid = function(createIntegrationGrid) {
    tt = createIntegrationGrid
    r = function(f, lowerLimit, upperLimit, ...) {
        n = length(lowerLimit)
        d = upperLimit - lowerLimit
        grid = createIntegrationGrid(n) # generates grid on [0, 1] ** n
        ## transform to [lowerLimit, upperLimit]
        nodes = sweep(sweep(grid[['nodes']], 2, d, '*'), 2, lowerLimit, '+')

        r = apply(nodes, 1, f, ...)
        r = (r %*% grid[['weights']])[1, 1] * prod(d)
        return(r)
    }
    return(r)
}

nint_integrateNCube = nint_integrateNCube_integrate(integrateA)


#' Integrate N Function
#'
#' Inferface to the integration over function dimensions.
#'
#' \code{nint_integrate} uses \code{nint_integrateNFunc} to handle function dimensions.
#' See examples below on how to deploy different solutions.
#'
#' @usage nint_integrateNFunc(f, funcs, x0, i0, ...)
#'
#' @param f the scalar-valued wrapper function to be integrated.
#' @param funcs the list of function dimensions.
#' @param x0 the partially realized point in the space.
#' @param i0 the vector of indices of function dimensions in the space.
#' @param ... other arguments passed to \code{f}.
#'
#' @return \code{nint_integrateNFunc} returns a single numeric.
#'
#' @seealso \code{\link{nint_integrate}}
#'
#' @example R/examples/nint_integrateNFunc.R
#'
#' @name nint_integrateNFunc
#'
#' @export
NULL

#' @details The function built by \code{nint_integrateNFunc_recursive} directly sums over discrete dimensions and uses \code{integrate1} otherwise.
#' In conjunction with \code{integrateA} this is the default.
#'
#' @param integrate1 \code{function(f, lowerLimit, upperLimit, ...)} which performs one dimensional integration.
#'
#' @return \code{nint_integrateNFunc_recursive} returns a recursive implementation for \code{nint_integrateNFunc}.
#'
#' @seealso \code{\link{integrateA}}
#'
#' @rdname nint_integrateNFunc
#'
#' @export
nint_integrateNFunc_recursive = function(integrate1) {
    r = function(f, funcs, x0, i0, ...) {
        maxDepth = length(funcs)
        if (maxDepth == 0)
            return(0)

        x = rep(0, maxDepth)
        d = 0
        g = Vectorize(function(xx, ...) {
            x[d] <<- xx
            x0[i0[d]] <<- xx
            if (d == maxDepth) {
                return(f(x, ...))
            }

            d <<- d + 1
            n = funcs[[d]](x0)
            type = nint_dtype(n)
            if (is.null(type) || type != nint_TYPE_INTV_DIM)
                r = sum(vapply(n, g, 0, ...))
            else
                r = integrate1(g, n[1], n[2], ...)
            d <<- d - 1
            return(r)
        }, 'xx')

        return(g(0, ...))
    }
    return(r)
}

nint_integrateNFunc = nint_integrateNFunc_recursive(function(...) integrateA(...)[['value']])



#' Integrate
#'
#' \code{nint_integrate} performs summation and integration of a scalar-valued function over a space or list structure of spaces.
#'
#' \code{nint_integrate} uses \code{nint_integrateNCube} and \code{nint_integrateNFunc} to handle interval and function dimensions.
#' See their help pages on how to deploy different solutions.
#'
#' The order of dimensions is optimized for efficiency.
#' Therefore interchangeability (except for function dimensions) is assumed.
#'
#' @param f the scalar-valued function (integrand) to be integrated.
#' @param space a space or list structure of spaces.
#' @param ... other arguments passed to \code{f}.
#'
#' @return \code{nint_integrate} returns a single numeric.
#'
#' @seealso \code{\link{nint_space}}, \code{\link{nint_transform}}, \code{\link{nint_integrateNCube}}, \code{\link{nint_integrateNFunc}}, \code{\link{fisherI}}
#'
#' @example R/examples/nint_integrate.R
#'
#' @export
nint_integrate = function(f, space, ...) {
    ispaces = nint_ispaces(space)
    if (length(ispaces) == 0)
        return(0)

    ## globals
    x = rep(NA_real_, sum(sapply(ispaces[[1]], function(x) length(x$i))))
    ispace = NULL
    gg = list() # sequence of functions
    maxDepth = 0
    d = 0 # depth
    i = 0 # indices
    g = NULL # data

    gs = function(...) sum(apply(g, 1, gd, ...)) # g scatter
    gi = function(...) nint_integrateNCube(gd, g[,1], g[,2], ...) # g interval
    gf = function(...) nint_integrateNFunc(gd, g, x, i, ...) # g function
    gl = list(s=gs, g=gs, i=gi, f=gf) # g list

    gd = function(xx, ...) { # g dispatch
        x[i] <<- xx
        if (d == maxDepth)
            return(f(x, ...))

        ## save state
        ti = i
        tg = g

        ## prepare descend
        d <<- d + 1
        idim = ispace[[d]]
        i <<- idim[['i']]
        g <<- idim[['g']]

        ## descend
        r = gg[[d]](...)

        ## restore state
        d <<- d - 1
        i <<- ti
        g <<- tg
        return(r)
    }

    af = f # argument
    r = 0
    for (ispace in ispaces) {
        gg = gl[names(ispace)]
        maxDepth = length(ispace)

        r = r + gd(0, ...)
    }

    return(r)
}



