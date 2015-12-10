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
#' @param contrasts   contrasts to use with the regression. See the \code{contrasts.arg}
#'                     of \code{model.matrix.default}
#' @param beta_init   initial beta vector to start at; when missing a first pass using
#'                     iolm is run to determine the starting point.
#' @param s_init      inital value of the scale factor; when missing a first pass using
#'                     iolm is run to determine the starting scale.
#' @param a           tuning parameter for Tukey bisquare. See Details for more information.
#' @param maxit       the limit on the number of IWLS iterations.
#' @param acc         accuracy for the IWLS stopping criterion.
#' @param trace       logical indicating if output should be produced for each
#'                     iteration.
#' @param tol         numeric tolerance. Set to -1 to ignore.
#' @details
#' The parameter \code{a} controls the tradeoff between efficency
#' and robustness. The default value of 4.685 yields an efficency of 95%
#' but breakdown point of about 10%; conversely 1.547 has only a 28% efficency
#' but a breakdown point of 50%. For large datasets, it may be preferable to set
#' the values of \code{a} lower than the default.
#'
#' @export
iorlm = function(formula, data, weights=NULL, subset=NULL,
                na.action=NULL, offset=NULL, contrasts=NULL,
                beta_init=NULL, s_init=NULL, a=4.685, maxit = 20,
                acc = 1e-4, trace=FALSE, tol=-1) {
  call <- match.call()

  if (!inherits(data, "adf")) data = adf(data)

  if (!is.null(weights) && !is.character(weights <- weights[[1]]))
    stop("weights must be a length one character vector")
  if (!is.null(subset) && !is.character(subset <- subset[[1]]))
    stop("subset must be a length one character vector")
  if (!is.null(offset) && !is.character(offset <- offset[[1]]))
    stop("offset must be a length one character vector")

  if (is.null(beta_init) || is.null(s_init))
    lm.out = iolm(formula, data, subset=subset, weights=weights,
                    na.action=na.action, offset=offset, contrasts=NULL,
                    tol=tol)

  beta <- if (is.null(beta_init)) as.numeric(lm.out$coefficients) else beta_init
  s <- if (is.null(s_init)) summary(lm.out)$sigma else s_init
  beta_old = beta

  converged=FALSE
  for (i in 1:maxit) {
    pvar = list(beta=beta, a=a, s=s)
    cvs = adf.apply(x=data, type="sparse.model",
                    FUN=rlm_kernel ,passedVars=pvar, formula=formula,
                    subset=subset,weights=weights, na.action=na.action,
                    offset=offset, contrasts=contrasts)
    cvs = cvs[!sapply(cvs,is.null)]

    XTWX = Reduce(`+`, Map(function(x) x$XTWX, cvs))
    XTWz = Reduce(`+`, Map(function(x) x$XTWz, cvs))
    resid20 = unlist(Map(function(x) x$resid20, cvs))

    s = median(abs(resid20 - median(resid20))) * 1.4826

    beta = Matrix::solve(XTWX, XTWz, tol=2*.Machine$double.eps)
    err = as.vector(Matrix::crossprod(beta-beta_old) / sum(beta_old^2))
    if (!is.null(beta_old) && err < acc) {
      converged=TRUE
      break
    }
    if (trace) cat(sprintf("Delta: %02.4f Scale: %02.4f Iterations - %d\n",err,s,i))
    beta_old = beta
  }

  b = as.numeric(beta)
  names(b) = rownames(beta)

  ret = list(coefficients=b,
    data=data,
    xtwx=XTWX,
    xtwz=XTWz,
    iter=i,
    a=a,
    s=s,
    converged=converged,
    formula=formula,
    call=call,
    num_obs=nobs)
  class(ret) = c("iorlm")
  ret
}

#' Print iorlm object
#'
#' @method print iorlm
#' @param x        output of iorlm
#' @param digits   significant digits to print
#' @param ...      optional, currently unused, arguments
#' @export
print.iorlm =
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

rlm_kernel = function(d, passedVars=NULL) {
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
  a <- if (!is.null(passedVars$a)) passedVars$a else 4.685
  s <- if (!is.null(passedVars$s)) passedVars$s else 1.0

  tbw = function(z, a) {
    out = (1 - (z/a)^2)^2
    out[abs(z) > a] = 0
    out
  }

  r = (d$x %*% passedVars$beta + offset - d$y)
  W = as.vector(weights * as.numeric(tbw(r / s, a = a)))

  list(XTWX=Matrix::crossprod(d$x, W * d$x),
       XTWz=Matrix::crossprod(d$x, W*d$y),
       cumulative_weight=sum(d$w),
       nobs=nobs,
       wy=Matrix::crossprod(sqrt(weights), d$y),
       resid20=quantile(r,seq(0,1,0.05)))
}
