#' Perform a linear regression
#'
#' @param formula   the formula for the regression
#' @param data      an abstract data frame, or something which can be
#'                  coerced to one.
#' @param subset    an options character string, which will be evaluated in the
#'                  frame of the data, to indicate which rows to include
#'                  in the analysis
#' @param weights   a optional character string, which will be evaluated in the
#'                  frame of the data, giving the sample weights for the regression
#' @param na.action a function which indicates what should happen when the data
#'                  contain 'NA's. See lm.fit for more details.
#' @param offset    a optional character string, which will be evaluated in the
#'                  frame of the data, giving the offsets for the regression
#' @param contrasts contrasts to use with the regression. See the \code{contrasts.arg}
#'                  of \code{model.matrix.default}
#' @param tol       numeric tolerance. Set to -1 to ignore.
#' @param parallel    integer. the number of parallel processes to use in the
#'                     calculation (*nix only).
#' @export
iolm = function(formula, data, subset=NULL, weights=NULL,
                na.action=NULL, offset=NULL, contrasts=NULL, tol=-1,
                parallel=1L) {
  call = match.call()
  if (!inherits(data, "adf")) data = as.adf(data)

  cvs = adf.apply(x=data, type="sparse.model",
          FUN=function(d,passedVars) {
            if (nrow(d$x) == 0L) return(NULL)
            z = d$y # non-offset version
            if (!is.null(d$offset)) d$y = d$y - d$offset
            if (!is.null(d$w)) {
              if (any(d$w == 0)) {
                ok = d$w != 0
                d$w = d$w[ok]
                d$x = d$x[ok,,drop = FALSE]
                d$y = d$y[ok]
                if (!is.null(d$offset)) d$offset = d$offset[ok]
              }
              sum_y = sum(d$y * d$w)
              sum_z = sum(z * d$w)
              sum_w = sum(d$w)
              if (!is.null(d$offset)) sum_o = sum(d$offset) else sum_o = 0
              d$x = d$x * sqrt(d$w)
              d$y = d$y * sqrt(d$w)
              if (!is.null(d$offset)) d$offset = d$offset * sqrt(d$w)
              z = z * sqrt(d$w)
            } else {
              sum_y = sum(d$y)
              sum_z = sum(z)
              sum_w = nrow(d$x)
            }
            if (!is.null(d$offset)) {
              oto = Matrix::crossprod(d$offset)
              xto = Matrix::crossprod(d$x, d$offset)
            } else {
              oto = Matrix(0)
              xto = Matrix(0,nrow=ncol(d$x))
            }
            return(list(xtx = Matrix::crossprod(d$x),
                        xty = Matrix::crossprod(d$x, d$y),
                        yty = Matrix::crossprod(d$y),
                        oto = oto,
                        xto = xto,
                        n = nrow(d$x),
                        sum_y = sum_y,
                        sum_w = sum_w,
                        sum_z = sum_z,
                        mean_x = apply(d$x,2,sum),
                        contrasts=attr(d$x, "contrasts"),
                        mt=d$mt))

          },formula=formula,subset=subset,weights=weights,
            na.action=na.action, offset=offset, contrasts=contrasts,
            parallel=parallel)

  cvs = cvs[!sapply(cvs,is.null)]
  if (length(cvs) == 0L) stop("No valid data.")
  xtx = Reduce(`+`, Map(function(x) x$xtx, cvs))
  xty = Reduce(`+`, Map(function(x) x$xty, cvs))
  sum_y = Reduce(`+`, Map(function(x) x$sum_y, cvs))
  sum_w = Reduce(`+`, Map(function(x) x$sum_w, cvs))
  yty = Reduce(`+`, Map(function(x) x$yty, cvs))
  oto = Reduce(`+`, Map(function(x) x$oto, cvs))
  xto = Reduce(`+`, Map(function(x) x$xto, cvs))
  sum_z = Reduce(`+`, Map(function(x) x$sum_z, cvs))
  n = Reduce(`+`, Map(function(x) x$n, cvs))
  mean_x = Reduce(`+`, Map(function(x) x$mean_x, cvs)) / n
  contrasts = cvs[[1]]$contrasts
  mt = cvs[[1]]$mt


  # Get rid of colinear variables.
  ch = chol(xtx, pivot=TRUE, tol=tol)
  if (attr(ch, "rank") < ncol(xtx)) {
    # xtx is rank deficient.
    new_name_order = attr(ch, "pivot")
    effective_rank = attr(ch, "rank")
    ft = mt
    formula_vars = colnames(attr(ft, "factors"))
    keep_vars = new_name_order[1:effective_rank]
    drop_var_names = colnames(xtx)[setdiff(new_name_order, keep_vars)]
    if (length(setdiff(drop_var_names, formula_vars)) > 0) {
      # A collinearity exists among factor levels in the design matrix.
      # Error out with an appropriate message.
      stop(paste("The following design matrix variables are colinear but",
                 "could not be dropped because they are associated with a",
                 "factor variable:",
                 setdiff(drop_var_names, formula_vars),
                 sep="\n"))
    } else {
      # The variable to drop can be dropped but we'll warn the user.
      warning(paste("Removing the following variables due to colinearity:",
              colnames(xtx)[tail(new_name_order,nrow(xtx)-effective_rank)],
              sep="\n"))
    }
    #Reorder xtx and xty in terms of the pivots up to the rank.
    xtx = xtx[keep_vars, keep_vars]
    xty = xty[keep_vars,]
    # Now update the formula so that it doesn't have the colinear terms.
    drop_inds = match(drop_var_names, formula_vars)
    ft = drop.terms(ft, drop_inds, keep.response=TRUE)
    form = formula(ft)
  }

  coefficients = as.numeric(Matrix::solve(xtx,xty))
  names(coefficients) = rownames(xtx)
  terms = mt
  contrasts = contrasts

  ret = list(coefficients=coefficients, call=call, terms=terms,
             xtx=xtx, xty=xty, yty=yty, oto=oto, xto=xto, sum_z=sum_z,
             mean_x=mean_x,
             sum_y=as.numeric(sum_y), sum_w=as.numeric(sum_w),
             n=n, data=data, contrasts=contrasts, rank=ncol(xtx))
  class(ret) = "iolm"
  ret
}

#' Get the regression diagnostics for a linear regression
#'
#' @method summary iolm
#' @param object   an object return from iolm
#' @param ...      optional, currently unused, arguments
#' @export
summary.iolm = function(object, ...) {
  call = match.call()

  beta = coef(object)
  btxtxb = t(beta) %*% object$xtx %*% beta
  rss = object$yty + btxtxb - 2 * t(beta) %*% object$xty

  mss = btxtxb + object$oto + 2 * Matrix::t(object$xto) %*% beta
  if (attr(object$terms, "intercept"))
    mss = mss - object$sum_z^2 / object$sum_w

  rss_beta = object$yty - btxtxb
  df = c(length(beta),
         object$n-length(beta),
         length(beta))
  sigma = as.numeric(sqrt(rss_beta / df[2]))
  #vcv = sigma * sqrt(solve(object$xtx))
  cov_unscaled = solve(object$xtx)
  se = sigma * sqrt(diag(cov_unscaled))
  tv = object$coefficients/se
  ptv = 2 * pt(abs(tv), df[2], lower.tail = FALSE)

  coefficients = cbind(object$coefficients, se, tv, ptv)
  colnames(coefficients) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(coefficients) = rownames(object$xtx)

  # Aliasing could be better.
  aliased = rep(FALSE, length(object$coefficients))
  names(aliased) = names(object$coefficients)
  df.int <- if (attr(object$terms, "intercept")) 1L else 0L
  r_squared = as.numeric(mss/(mss + rss))
  adj_r_squared = 1 - (1-r_squared) * sum(df[1:2], -1*df.int)/df[2]
  fstatistic = c((mss/(df[1]-df.int))/sigma^2, df[1]-df.int, dendf=df[2])
  names(fstatistic) = c("value", "numdf", "dendf")

  ret = list(call=call, terms=object$terms, coefficients=coefficients,
       aliased=aliased, sigma=sigma, df=df, r.squared=r_squared,
       adj.r.squared=as.numeric(adj_r_squared),
       fstatistic=sapply(fstatistic,as.numeric),
       cov.unscaled=cov_unscaled, data=object$data, dfpp=object$dfpp)
  class(ret) = "summary.iolm"
  ret
}

#' Print iolm object
#'
#' @method print iolm
#' @param x        output of iolm
#' @param digits   significant digits to print
#' @param ...      optional, currently unused, arguments
#' @export
print.iolm =
function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L,
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

#' Print iolm summary
#'
#' @method print summary.iolm
#' @param x             output of summary.iolm
#' @param digits        significant digits to print
#' @param symbolic.cor  logical. Should symbolic correlation be printed.
#' @param signif.stars  logical. Should significant starts be printed.
#' @param ...           optional, currently unused, arguments
#' @export
print.summary.iolm =
function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  if (length(x$aliased) == 0L) {
      cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L])
        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
            sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
            colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
        na.print = "NA", ...)
  }
  cat("\nResidual standard error:", format(signif(x$sigma,
      digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  if (nzchar(mess <- naprint(x$na.action)))
      cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared,
        digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L],
        digits = digits), "on", x$fstatistic[2L], "and",
        x$fstatistic[3L], "DF,  p-value:", format.pval(pf(x$fstatistic[1L],
            x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE),
            digits = digits))
    cat("\n")
  }
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2,
          digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}
