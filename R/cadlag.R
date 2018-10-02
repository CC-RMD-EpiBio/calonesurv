#  Derived from File src/library/stats/R/stepfun.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2014 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## Constructor for  Step Functions:

## Given x[1L] .. x[n] ;	 y[0] .. y[n]  (one value more !)
## For 'cadlag' functions :  f(t) = y[i]  iff  t in  ( x[i], x[i+1] ],
##     where  x[0] := - Inf
##
## 'General case':  f(x[i]) = z[i]  with arbitrary z[]
## -- but we would have to modify 'approxfun' or not be able to use it..
## --->> forget about general case

#' @export
setClass("cadlag",contains ="function",representation(order = "numeric", x = "numeric", y = "numeric"))

#' Cadlag linear functions
#'
#' This function allows you to express your love of cats.
#' @param x c() of x values
#' @param y  c() of y values of length 1 greater than x
#' @keywords cadlag
#' @export
#' @examples
#' cadlag(c(0,1),c(1,2,3))
#'
cadlag = function(x,y=FALSE, order = 0){

  if(class(x) == "survfit"){
    t = x$time
    y = x$surv
    if(sum(abs(diff(y)))==0) return(1)
    obj = new("cadlag", cadlag.fun(t,c(1,y)), order = order, x = t, y = c(1,y))
  }

  else{
    dy = diff(y)
    if(sum(na.omit(dy)==0)>0){
      xnew = x[dy!=0]
      ynew = y[-(which(dy==0) + 1L)]

      if(length(xnew)==0){
        xnew = c(x[1],max(x))
        ynew = y[1:3]
      }
    }
    else{
      xnew = x
      ynew = y
    }
    #if(sum(abs(na.omit(diff(ynew))))==0) return(ynew[1])
    obj = new("cadlag", cadlag.fun(xnew,ynew, order = order), order = order, x = xnew, y = ynew)
  }

  return(obj)
}

#' Cadlag linear functions
#'
#' This function allows you to express your love of cats.
#' @param x c() of x values
#' @param y  c() of y values of length 1 greater than x
#' @param order 0 for constant and 1 for linear
#' @keywords cadlag
#' @export
#' @examples
#' cadlag.fun(c(0,1),c(1,2,3))
#' cadlag.fun(c(0,1),c(1,2),order=1)
#'
cadlag.fun = function(x, y, order = 0)
{
  if(is.unsorted(x)) stop("cadlag: 'x' must be ordered increasingly")
  if(length(x)!=length(unique(x))) stop("cadlag: 'x' must consist of unique entries")
  n <- length(x)
  if(n < 1) stop("'x' must have length >= 1")
  if(order == 0){
    n1 <- n + 1L
    if(length(y) != n1) stop("'y' must be one longer than 'x'")
    rval = approxfun(x, y[-1], method = "constant",
                     yleft = y[1L], yright = y[n1], f = 0, ties = "ordered")
  }
  else if(order == 1){
    if(length(y) != length(x)) stop("'y' must the same size as 'x'")
    rval = approxfun(x, y, method = "linear", rule = 1,
                     f = 0, ties = "ordered")
  }
  else{
    stop("Order must be 0 or 1")
  }
  rval
}

is.cadlag <- function(x) inherits(x, "cadlag")

as.cadlag <- function(x, ...) UseMethod("as.cadlag")
as.cadlag.default <- function(x, ...)
{
  if(is.cadlag(x)) x
  else stop("no 'as.cadlag' method available for 'x'")
}

## Quite obvious  that I will want to have  knots.spline(..)  etc......
knots         <- function(Fn, ...) UseMethod("knots")
knots.cadlag <- function(Fn, ...) eval(expression(x), envir=environment(Fn))


print.cadlag <- function (x, digits = getOption("digits") - 2, ...)
{
  numform <- function(x) paste(formatC(x, digits = digits), collapse=", ")
  i1 <- function(n) 1L:min(3L, n)
  i2 <- function(n) if(n >= 4L) max(4L, n-1L):n else integer()
  cat("Cadlag function: ")
  env <- environment(x)
  n <- length(xx <- eval(expression(x), envir = env))
  cat(" x[1:", n, "] = ", numform(xx[i1(n)]),
      if(n > 3L) ", ", if(n > 5L) " ..., ", numform(xx[i2(n)]), "\n", sep = "")
  y <- eval(expression(c(yleft, y)), envir = env)
  cat(n+1L, " plateau levels = ", numform(y[i1(n+1L)]),
      if(n+1L > 3L) ", ", if(n+1L > 5L) " ..., ", numform(y[i2(n+1L)]), "\n",
      sep = "")
  invisible(x)
}

summary.cadlag <- function(object, ...)
{
  n <- length(eval(expression(x), envir = environment(object)))
  if(!is.integer(n) || n < 1L) stop("not a valid step function")
  cat("Step function with continuity 'f'=",
      format(eval(expression(f), envir = environment(object))),
      ", ", n, if(n <= 6L) "knots at\n" else "knots with summary\n")
  summ <- if(n > 6L) summary else function(x) x
  print(summ(knots(object)))
  cat(if(n > 6L) "\n" else "  ", "and	", n+1L,
      " plateau levels (y) ", if(n <= 6L) "at\n" else "with summary\n",
      sep  = "")
  print(summ(eval(expression(c(yleft,y)), envir = environment(object))))
  invisible()
}

## Purpose: plot method for  cadlag (step function) objects
## --------------------------------------------------------------------
## Arguments: for numeric 'x', do empirical CDF;	  ==> `` ?plot.step ''
## --------------------------------------------------------------------
## Author: Martin Maechler <maechler@stat.math.ethz.ch>
##	      1990, U.Washington, Seattle; improved, Dec.1993
##	      Ported to R :  Sept.1997.
plot.cadlag <-
  function(x, xval, xlim, ylim = range(c(y,Fn.kn)),
           xlab = "x", ylab = "f(x)", main = NULL,
           add = FALSE, verticals = TRUE, do.points = (n < 1000),
           pch = par("pch"), col = par("col"),
           col.points = col, cex.points = par("cex"),
           col.hor = col, col.vert = col,
           lty = par("lty"), lwd = par("lwd"),
           ...)
  {
    if(!is.cadlag(x)) { #- make it work when called explicitly with data
      if(is.numeric(x)) {
        sarg <- substitute(x)
        x <- ecdf(x)
        attr(x,"call") <- call("ecdf", sarg)
      } else stop("'plot.cadlag' called with wrong type of argument 'x'")
    }
    order = slot(x,"order")
    if(missing(main))
      main <- {
        cl <- attr(x,"call")
        deparse(if(!is.null(cl))cl else sys.call())
      }

    knF <- knots(x)
    xval <- if (missing(xval)) knF else sort(xval)
    if (missing(xlim)) {
      rx <- range(xval)
      dr <-
        if(length(xval) > 1L)
          max(0.08 * diff(rx), median(diff(xval)))
      else
        abs(xval)/16
      xlim <- c(slot(x,'x')[1],slot(x,'x')[length(slot(x,'x'))])

    } else dr <- diff(xlim)

    xval <- xval[xlim[1L]-dr <= xval & xval <= xlim[2L]+dr]

    ## Careful for heights of horizontals -- these depend on f
    ti <- c(xlim[1L]-dr, xval, xlim[2L]+dr)
    ti.l <- ti[-length(ti)]
    ti.r <- ti[-1L]
    y <- x(0.5*(ti.l + ti.r))
    n <- length(y)
    Fn.kn <- x(xval)

    ##------------------------ Plotting ----------------------------

    dev.hold(); on.exit(dev.flush())
    ## horizontal segments
    if (add)
      segments(ti.l, y, ti.r, y, col=col.hor, lty=lty, lwd=lwd, ...)
    else {
      if(missing(ylim)) ylim <- range(na.omit(c(y,Fn.kn)))
      #print(na.omit(Fn.kn))
      #print(y)
      plot(NA, NA, type = "n", xlim = xlim, ylim = ylim,
           xlab = xlab, ylab = ylab, main = main, ...)
      segments(ti.l, y, ti.r, y, col = col.hor, lty = lty, lwd = lwd)
    }
    if(do.points){ #points(xval, Fn.kn, pch = pch, col = col.points, cex = cex.points)
      Fn.kn.1 = slot(x,'y')
      points(xval, Fn.kn.1[-length(Fn.kn.1)], pch = 1, col = col.points, cex = cex.points)
      points(xval, Fn.kn, pch = 19, col = col.points, cex = cex.points)
    }


    if(verticals)
      segments(xval, y[-n], xval, y[-1L], col = col.vert, lty = lty, lwd = lwd)
    invisible(list(t = ti, y = y))
  }

lines.cadlag <- function(x, ...) plot(x, add = TRUE, ...)

as.cadlag.isoreg <- function(x, ...)
{
  sf <- cadlag(x = (if(x$isOrd) x$x else x$x[x$ord])[x$iKnots],
               y = c(x$yf[x$iKnots], x$yf[length(x$yf)]))
  attr(sf, "call") <- x$call
  sf
}


setMethod(Math, "cadlag", function(x) {
  y = slot(x,'y')
  ynew = getFunction(.Generic)(y)
  return(cadlag(slot(x,'x'),ynew))
})

#setMethod(Ops, "cadlag", function(e1, e2){
# default OPS pass through
#  return(getFunction(.Generic)(slot(e1,'y'),slot(e2,'y')))
#})

setMethod("+", signature(e1 = "cadlag", e2 = "numeric"), function(e1, e2) {
  y = slot(e1,'y')+e2
  x = slot(e1,'x')
  return(cadlag(x,y,order = slot(e1,'order')))
})

setMethod("+", signature(e1 = "numeric", e2 = "cadlag"), function(e1, e2) {
  return(e2+e1)
})

setMethod("-", signature(e1 = "numeric", e2 = "cadlag"), function(e1, e2) {
  return(e1+ -1*e2)
})

setMethod("*", signature(e1 = "numeric", e2 = "cadlag"), function(e1, e2) {
  y = slot(e2,'y')*e1
  x = slot(e2,'x')
  return(cadlag(x,y,order = slot(e2,'order')))
})

setMethod("*", signature(e1 = "cadlag", e2 = "numeric"), function(e1, e2) {
  return(e2*e1)
})

setMethod("/", signature(e1 = "numeric", e2 = "cadlag"), function(e1, e2) {
  if(slot(e2,'order')>0) stop("not implemented")
  y = e1/slot(e2,'y')
  x = slot(e2,'x')
  return(cadlag(x,y))
})

setMethod("/", signature(e1 = "cadlag", e2 = "numeric"), function(e1, e2) {
  return(e1*(1/e2))
})

setMethod("/", signature(e1 = "cadlag", e2 = "cadlag"), function(e1, e2) {
  if(slot(e2,'order')>0) stop("not implemented")
  return(e1*(1/e2))
})

setMethod("+", signature(e1 = "cadlag", e2 = "cadlag"), function(e1, e2) {
  t1 = slot(e1,'x')
  t2 = slot(e2,'x')
  t = sort(unique(c(t1,t2)),decreasing = F)
  y = e1(t) + e2(t)

  order1 = slot(e1,"order")
  order2 = slot(e2,"order")
  order = max(order1,order2)

  if(order ==0)
    y = c(e1(t[1]-1000)+e2(t[1]-1000),y)

  return(cadlag(t,y,order=order))
})

setMethod("-", signature(e1 = "cadlag", e2 = "cadlag"), function(e1, e2) {
  return(e1 + -1*e2)
})

setMethod("-", signature(e1 = "cadlag", e2 = "numeric"), function(e1, e2) {
  return(e1 + -e2)
})

setMethod("*", signature(e1 = "cadlag", e2 = "cadlag"), function(e1, e2) {
  order1 = slot(e1,"order")
  order2 = slot(e2,"order")
  order = max(order1,order2)
  t1 = slot(e1,'x')
  y1 = slot(e1,'y')
  t2 = slot(e2,'x')
  y2 = slot(e2,'y')
  t = sort(unique(c(t1,t2)),decreasing = F)

  if(order == 0){
    y = e1(t) * e2(t)
    y = c(y1[1]*y2[1],y)
  }
  else if(order1 + order2 == 1){
    y = e1(t) * e2(t)
  }
  else{
    stop("Not implemented yet")
  }
  return(cadlag(t,y,order = order))
})


setMethod("==", signature(e1 = "cadlag", e2 = "numeric"), function(e1, e2) {
  y = slot(e1,'y')
  return(sum(y != e2) == 0)
})

setMethod("==", signature(e1 = "numeric", e2 = "cadlag"), function(e1, e2) {
  return(e2 == e1)
})


# definite and indefinite integral

#environment(integrate.cadlag) = environment(cadlag)

#' Indefinite integral of cadlag functions
#'
#' This function allows you to express your love of cats.
#' @param f object of type cadlag
#' @return object of type cadlag
#' @keywords cadlag
#' @export
#' @examples
#' f = cadlag(c(0,1),c(1,2,3))
#' F = integrate.cadlag(f)
#'
integrate.cadlag = function (f, lower = NA, upper = NA, ..., subdivisions = 100L, rel.tol = .Machine$double.eps^0.25,
                             abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE,
                             aux = NULL)
{
  y = slot(f,"y")
  x = slot(f,"x")
  dy = diff(y)
  dx = diff(x)

  if(is.na(lower) | is.na(upper)){
    if(slot(f,"order")>0) stop("Indefinite integration only implemented for peicewise constant functions")

    # we will return a function with the indefinite integral
    masses = dx * y[-c(1,length(y))]
    ynew = c(0,cumsum(masses))
    return(cadlag(x,ynew,order=1))
  }
  else{

  }

}

#setMethod("integrate",signature(f = "cadlag"),integrate.cadlag)

#setGeneric("integrate",integrate.cadlag, signature="cadlag")

### Testing

test_cadlag = new.env()
with(test_cadlag,{
  f1 = cadlag(1:10,sample(11))

  x = 1:6
  y = c(1,1,1,2,3,4,5)
  f2 = cadlag(x,y)

})

integrate.cadlag_measure = function(f,F, lower, upper){
  # F is a cadlag function of order 0
  if(slot(F,'order')!=0) stop("F must be peicewise constant")
  x_F = slot(F,'x')
  x_F = x_F[x_F <= upper & x_F >=lower]
  eps = 1e-12
  dF = F(x_F) - F(x_F-eps)
  return(sum(f(x_F)*dF))
}
