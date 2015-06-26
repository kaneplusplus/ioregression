#' Fit an Lp Regression
#'
#' Minimizes the loss ||Y - Xb||_p, for a suitable choice
#' of p.
#'
#' @importFrom       Matrix crossprod
#' @param formula     the formula for the regression
#' @param data        an abstract data frame, or something which can be
#'                     coerced to one.
#' @param p           value of p for the regression, such that 1 <= p < 2.
#'                     Using p=2 is also allowed for testing purposes, but
#'                     sub-optimal as iolm should instead be called directly.
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
#' @param beta_init   initial beta vector to start at; when missing an first pass using
#'                     iolm is run to determine the starting point.
#' @param delta       tuning parameter to provide weights for very small residuals.
#'                     most users should not need to change this parameter.
#' @param maxit       the limit on the number of IWLS iterations.
#' @param acc         accuracy for the IWLS stopping criterion.
#' @param trace       logical indicating if output should be produced for each
#'                     iteration.
#' @param tol         numeric tolerance. Set to -1 to ignore.
#' @export
iolp = function(formula, data, p=1.5, weights=NULL, subset=NULL,
                na.action=NULL, offset=NULL, contrasts=NULL,
                beta_init=NULL, delta=0.0001, maxit=20, acc=1e-5,
                trace=FALSE, tol=-1) {
  call <- match.call()

  if (p > 2 || p < 1) stop("Invalid input to p; must have 1 <= p <= 2.")
  if (!inherits(data, "adf")) data = as.adf(data)

  if (is.null(beta_init)) {
    lm.out = iolm(formula, data, subset=subset, weights=weights,
                    na.action=na.action, offset=offset, contrasts=NULL,
                    tol=tol)
    beta = beta_old = as.numeric(lm.out$coefficients)
  } else beta = beta_old = beta_init

  converged=FALSE
  for (i in 1:maxit) {
    pvar = list(beta=beta, p=p, delta=delta)
    cvs = adf.apply(x=data, type="sparse.model",
                    FUN=lp_kernel, passedVars=pvar, formula=formula,
                    subset=subset, weights=weights, na.action=na.action,
                    offset=offset, contrasts=contrasts)
    cvs = cvs[!sapply(cvs,is.null)]


    XTWX = Reduce(`+`, Map(function(x) x$XTWX, cvs))
    XTWz = Reduce(`+`, Map(function(x) x$XTWz, cvs))

    beta = Matrix::solve(XTWX, XTWz, tol=2*.Machine$double.eps)
    err = as.vector(Matrix::crossprod(beta-beta_old) / sum(beta_old^2))
    if (!is.null(beta_old) && err < acc) {
      converged=TRUE
      break
    }
    if (trace) cat(sprintf("Delta: %02.4f Iterations - %d\n",err,i))
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
