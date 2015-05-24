#' Fit a robust linear regression
#'
#' Fits an M-estimator using Tukey's bisquare function for
#' estimating slope coefficients in a linear model.
#'
#' @param formula     the formula for the regression
#' @param data        an abstract data frame, or something which can be
#'                     coerced to one.
#' @param weights     a optional character string, which will be evaluated in the
#'                     frame of the data, giving the sample weights for the regression
#' @param subset      an options character string, which will be evaluated in the
#'                     frame of the data, to indicate which rows to include
#'                     in the analysis
#' @param na.action   a function which indicates what should happen when the data
#'                     contain 'NA's. See lm.fit for more details.
#' @param offset      a optional character string, which will be evaluated in the
#'                     frame of the data, giving the offsets for the regression
#' @param a           tuning parameter for Tukey bisquare. See Details for more information.
#' @param maxit       the limit on the number of IWLS iterations.
#' @param acc         accuracy for the IWLS stopping criterion.
#' @param contrasts   contrasts to use with the regression. See the \code{contrasts.arg}
#'                     of \code{model.matrix.default}
#' @param tol         numeric tolerance. Set to -1 to ignore.
#' @details
#' The parameter \code{a} controls the tradeoff between efficency
#' and robustness. The default value of 4.685 yields an efficency of 95%
#' but breakdown point of about 10%; conversely 1.547 has only a 28% efficency
#' but a breakdown point of 50%. For large datasets, it may be preferable to set
#' the values of \code{a} lower than the default.
#'
#' @export
iolp = function(formula, data, p=1.5, weights=NULL, subset=NULL,
                na.action=NULL, offset=NULL, contrasts=NULL,
                beta_init=NULL, delta=0.0001, maxit=20, acc=1e-5, tol=-1) {
  call <- match.call()

  if (!inherits(data, "adf")) data = as.adf(data)

  if (is.null(beta_init)) {
    lm.out = iolm(formula, data, subset=subset, weights=weights,
                    na.action=na.action, offset=offset, contrasts=NULL,
                    tol=tol)
    beta = beta_old = as.numeric(lm.out$coefficients)
  } else beta = beta_old = beta_init

  converged=FALSE
  for (i in 1:maxit) {
    print(as.numeric(beta))

    pvar = list(beta=beta, p=p, delta=delta)
    cvs = adf.apply(x=data, type="sparse.model",
                    FUN=lp_kernel, passedVars=pvar, formula=formula,
                    subset=subset, weights=weights, na.action=na.action,
                    offset=offset, contrasts=contrasts)
    cvs = cvs[!sapply(cvs,is.null)]


    XTWX = Reduce(`+`, Map(function(x) x$XTWX, cvs))
    XTWz = Reduce(`+`, Map(function(x) x$XTWz, cvs))

    beta = Matrix::solve(XTWX, XTWz, tol=2*.Machine$double.eps)
    if (!is.null(beta_old) &&
        as.vector(Matrix::crossprod(beta-beta_old) / sum(beta_old^2 + delta)) < acc) {
      converged=TRUE
      break
    }
    beta_old = beta
  }

  b = as.numeric(beta)
  names(b) = rownames(beta)

  ret = list(coefficients=b,
    data=data,
    xtwx=XTWX,
    xtwz=XTWz,
    iter=i,
    p=p,
    converged=converged,
    formula=formula,
    call=call,
    num_obs=nobs)
  class(ret) = c("iolp")
  ret
}

#' Print iolp object
#'
#' @method print iolp
#' @param x        output of iolp
#' @param digits   significant digits to print
#' @param ...      optional, currently unused, arguments
#' @export
print.iolp =
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

lp_kernel = function(d, passedVars=NULL) {
  if (nrow(d$x) == 0L) return(NULL)
  if (!is.null(d$w)) {
    if (any(d$w == 0)) {
      ok = d$w != 0
      d$w = d$w[ok]
      d$x = d$x[ok,,drop = FALSE]
      d$y = d$y[ok]
      if (!is.null(d$offset)) d$offset = d$offset[ok]
    }
  }
  nobs = length(d$y)
  offset <- if (!is.null(d$offset)) d$offset else offset = rep.int(0,nobs)
  weights <- if (!is.null(d$w)) d$w else rep(1, nobs)
  p <- if (!is.null(passedVars$a)) passedVars$a else 1.5
  delta <- if (!is.null(passedVars$s)) passedVars$s else 0.0001

  r = (d$x %*% passedVars$beta + offset - d$y)
  Winv = abs(r)^(2-passedVars$p)
  Winv[Winv < passedVars$delta] = passedVars$delta
  W = as.vector(weights * 1/Winv)

  list(XTWX=Matrix::crossprod(d$x, W*d$x),
       XTWz=Matrix::crossprod(d$x, W*d$y),
       cumulative_weight=sum(d$w),
       nobs=nobs,
       wy=Matrix::crossprod(sqrt(weights), d$y))
}
