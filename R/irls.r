
ioirls = function(formula, family, data, weights, subset,
                na.action, start, etastart, mustart, offset,
                control, contrasts, trace, tol, beta_update_fun) {
  call <- match.call()
  control <- do.call("glm.control", control)
  if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
      family <- family()
  if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
  }

  if (!is.null(weights) && !is.character(weights <- weights[[1]]))
    stop("weights must be a length one character vector")
  if (!is.null(subset) && !is.character(subset <- subset[[1]]))
    stop("subset must be a length one character vector")
  if (!is.null(offset) && !is.character(offset <- offset[[1]]))
    stop("offset must be a length one character vector")

  converged=FALSE
  beta = beta_old = start
  wtdmu = deviance =  cumulative_weight = NULL
  for (i in 1:control$maxit) {
    pvar = list(beta=beta, cw=cumulative_weight, family=family,
                deviance=deviance, wtdmu=wtdmu)
    cvs = adf.apply(x=data, type="sparse.model", FUN=glm_kernel,
                    args=pvar, formula=formula, subset=subset,
                    weights=weights, na.action=na.action,
                    offset=offset, contrasts=contrasts)
    cvs = cvs[!sapply(cvs, is.null)]
    num_rows = Reduce(`+`, Map(function(x) x$num_rows, cvs))
    XTWX = Reduce(`+`, Map(function(x) x$XTWX, cvs))
    XTWz = Reduce(`+`, Map(function(x) x$XTWz, cvs))
    deviance = Reduce(`+`, Map(function(x) x$deviance, cvs))
    cumulative_weight = Reduce(`+`, Map(function(x) x$cumulative_weight, cvs))
    wtdmu = if(attributes(terms(formula))$intercept) {
      Reduce(`+`, Map(function(x) x$wy, cvs))/cumulative_weight
    } else {
      family$linkinv(0)
    }
    contrasts = cvs[[1]]$contrasts
    wx_norm = Reduce(`+`, Map(function(x) x$wx_norm, cvs))
    # TODO: Add checking for singularities here.
    beta = do.call("beta_update_fun")
    if (!is.null(beta_old))
      err = as.vector(Matrix::crossprod(beta-beta_old))
    else
      err = control$epsilon * 2
    p_val=Inf
    if ( (!is.null(beta_old) && (err < control$epsilon)) ) {
      converged=TRUE
      break
    }
    if (trace) {
      cat(sprintf("Delta: %02.4f Deviance: %02.4f Iterations - %d\n",
                  err,deviance,i))
    }
    beta_old = beta
  }

  # We only calculate these here, as we only care about the
  #  converging loop:
  aic = Reduce(`+`, Map(function(x) x$aic, cvs))
  RSS = Reduce(`+`, Map(function(x) x$RSS, cvs))
  nobs = Reduce(`+`, Map(function(x) x$nobs, cvs))
  null_dev = Reduce(`+`, Map(function(x) x$null_dev, cvs))
  beta_names = row.names(beta)
  beta = as.vector(beta)
  names(beta) = beta_names
  rank=nrow(XTWX)

  aic = aic + 2 * nrow(XTWX)

  # This will have to change when we support weights.
  nulldf = nobs - attributes(terms(formula))$intercept
  resdf = nobs - rank

  var_res = RSS/resdf
  dispersion = if (family$family %in% c("poisson", "binomial")) 1 else var_res

  ret = list(coefficients=beta,
    family=family,
    deviance = deviance,
    aic=aic,
    data=data,
    rank=rank,
    xtwx=XTWX,
    xtwz=XTWz,
    iter=i,
    dispersion=dispersion,
    rss=RSS,
    converged=converged,
    formula=formula,
    call=call,
    num_obs=nobs,
    df.null=nulldf,
    null.deviance=null_dev,
    df.residual=resdf,
    terms=terms.formula(formula),
    control=control,
    contrasts=contrasts)
  class(ret) = c("ioglm", "iolm")
  ret
}

glm_kernel = function(d, passedVars=NULL) {
  if (nrow(d$x) == 0L) return(NULL)
  if (!is.null(d$w)) {
    if (any(d$w == 0)) {
      ok = d$w != 0
      d$w = d$w[ok]
      d$x = d$x[ok,,drop = FALSE]
      d$y = d$y[ok]
      if (!is.null(d$offset)) d$offset = d$offset[ok]
    }
    sum_y = sum(d$y * d$w)
    sum_w = sum(d$w)
    # d$x = d$x * sqrt(d$w)
    # d$y = d$y * sqrt(d$w)
  } else {
    sum_y = sum(d$y)
    sum_w = nrow(d$x)
  }
  nobs = length(d$y)
  offset <- if (!is.null(d$offset)) d$offset else offset = rep.int(0,nobs)
  weights <- if (!is.null(d$w)) d$w else rep(1, nobs)
  family = passedVars$family

  if (is.null(passedVars$beta)) {
    # It's the first iteration.
    etastart = start = mustart = NULL
    y = d$y
    eval(family$initialize)
    eta = family$linkfun(mustart)
  } else {
    eta = as.numeric(d$x %*% passedVars$beta)
  }
  g = family$linkinv(eta <- eta + offset)
  gprime = family$mu.eta(eta)
  residuals = (d$y-g)/gprime
  z = eta - offset + residuals
  W = as.vector(weights * gprime^2 / family$variance(g))
  aic = NULL
  null_dev = NULL
  if (!is.null(passedVars$cw) && !is.null(passedVars$wtdmu)) {
    null_dev = sum(family$dev.resids(d$y, passedVars$wtdmu, weights))
  }
  RSS = sum(W*residuals^2)
  deviance = sum(family$dev.resids(d$y, g, weights))
  if (!is.null(passedVars$deviance) && !is.null(passedVars$cw)) {
    aic = family$aic(d$y, length(d$y), g, weights, deviance)
  }

  list(XTWX=Matrix::crossprod(d$x, W * d$x), XTWz=Matrix::crossprod(d$x, W*z),
       deviance=deviance, null_dev=null_dev, cumulative_weight=sum(d$w), 
       nobs=nobs, aic=aic, RSS=RSS, contrasts=attr(d$x, "contrasts"),
       wy=Matrix::crossprod(sqrt(weights), d$y), num_row=nrow(d$x))
}

#' Print ioglm object
#'
#' @method print ioglm
#' @param x        output of iolm
#' @param ...      optional arguments passed to print.glm
#' @export
print.ioglm =
function (x, ...) {
  class(x) <- "glm"
  print(x)
}

#' Get the regression diagnostics for a linear regression
#'
#' @method summary ioglm
#'
#' @param object    an object return from ioglm
#' @param ...       optional, currently unused, arguments
#' @export
summary.ioglm = function(object, ...) {
  call = object$call
  terms = object$terms
  dispersion = object$dispersion
  inv_scatter = solve(object$xtwx)
  standard_errors = sqrt(dispersion * diag(inv_scatter))
  stat_vals = object$coefficients/standard_errors
  if (object$family$family %in% c("binomial", "poisson")) {
    p_vals = 2 * pnorm(abs(stat_vals), lower.tail=FALSE)
  } else {
    p_vals = 2 * pt(abs(stat_vals), df=object$df, lower.tail=FALSE)
  }
  coef_names = names(object$coefficients)
  coefficients = cbind(object$coefficients, standard_errors, stat_vals, p_vals)
  colnames(coefficients) = c("Estimate", "Std. Error", "t value", "Pr(>|z|)")
  rownames(coefficients) = coef_names
  if (object$family$family %in% c("binomial", "poisson"))
    colnames(coefficients)[3] = "z value"
  aliased = object$coefficients
  aliased = FALSE
  ret = list(call=call,
       terms=terms,
       family=object$family,
       deviance=object$deviance,
       aic=object$aic,
       contrasts=object$contrasts,
       df.residual=object$df.residual,
       null.deviance=object$null.deviance,
       df.null=object$df.null,
       iter=object$iter,
       deviance.resid=NA,
       coefficients=coefficients,
       aliased=aliased,
       dispersion=object$dispersion,
       df=c(ncol(object$xtwx), object$df.residual, ncol(object$xtwx)),
       data=object$data,
       cov.unscaled=inv_scatter,
       cov.scaled=inv_scatter * dispersion)
  class(ret) = c("summary.glm")
  ret
}
