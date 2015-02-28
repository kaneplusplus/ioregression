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
#' @export
iolm = function(formula, data, subset=NULL, weights=NULL,
                na.action=NULL, offset=NULL, contrasts=NULL, tol=-1) {
  call = match.call()
  if (!inherits(data, "adf")) data = as.adf(data)

  cvs = adf.apply(x=data, type="sparse.model",
          FUN=function(d) {
            if (nrow(d$x) == 0L) return(NULL)
            return(list(xtx = Matrix::crossprod(d$x), xty = Matrix::crossprod(d$x, d$y),
                        xto = if (!is.null(offset)) Matrix::crossprod(d$x, d$y-d$offset) else NULL,
                        yty = Matrix::crossprod(d$y), n = nrow(d$x), sum_y = sum(d$y)))
          },formula=formula,subset=subset,weights=weights,
            na.action=na.action, offset=offset, contrasts=contrasts)
  cvs = cvs[!sapply(cvs,is.null)]
  if (length(cvs) == 0L) stop("No valid data.")
  xtx = Reduce(`+`, Map(function(x) x$xtx, cvs))
  xty = Reduce(`+`, Map(function(x) x$xty, cvs))
  xto = if (!is.null(offset)) Reduce(`+`, Map(function(x) x$xto, cvs)) else 0.0
  sum_y = Reduce(`+`, Map(function(x) x$sum_y, cvs))
  yty = Reduce(`+`, Map(function(x) x$yty, cvs))
  n = Reduce(`+`, Map(function(x) x$n, cvs))
  contrasts = cvs[[1]]$contrasts

  design_matrix_names = colnames(xtx)

  # Get rid of colinear variables.
  ch = chol(xtx, pivot=TRUE, tol=tol)
  if (attr(ch, "rank") < ncol(xtx)) {
    # xtx is rank deficient.
    new_name_order = attr(ch, "pivot")
    effective_rank = attr(ch, "rank")
    ft = terms(formula)
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
  coefficients = as.numeric(Matrix::solve(xtx,xty-xto)) #as.numeric(solve(xtx) %*% xty)
  names(coefficients) = rownames(xtx)
  terms = terms(formula)
  contrasts = contrasts

  ret = list(coefficients=coefficients, call=call, terms=terms,
             xtx=xtx, xty=xty, xto=xto, yty=yty,
             sum_y=as.numeric(sum_y), n=n, data=data,
             contrasts=contrasts, rank=ncol(xtx))
  class(ret) = "iolm"
  ret
}

#' Get the regression diagnostics for a linear regression
#'
#' @param object an object return from iolm
#' @export
summary.iolm = function(object, ...) {
  call = match.call()

  beta = coef(object)
  btxtxb = t(beta) %*% object$xtx %*% beta
  rss = object$yty + btxtxb - 2 * t(beta) %*% object$xty
  tss = object$yty - object$sum_y^2 / object$n
  mss <- if (attr(object$terms, "intercept"))
            (tss - rss) else as.numeric(btxtxb)

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
#' @param x        output of iolm
#' @param digits   significant digits to print
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
#' @param x             output of summary.iolm
#' @param digits        significant digits to print
#' @param symbolic.cor  logical. Should symbolic correlation be printed.
#' @param signif.stars  logical. Should significant starts be printed.
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
