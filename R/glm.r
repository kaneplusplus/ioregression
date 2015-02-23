
glm_kernel = function(y, mm, family, beta, deviance,
                      cumulative_weight, total_obs) {
  nobs = length(y)
  weights = rep(1, nobs)
  if (is.null(beta)) {
    # It's the first iteration.
    etastart = start = mustart = NULL
    eval(family$initialize)
    eta = family$linkfun(mustart)
  } else {
    eta = mm %*% beta
  }
  g = family$linkinv(eta)
  gprime = family$mu.eta(eta)
  residuals = (y-g)/gprime
  z = eta + residuals
  W = as.vector(gprime^2 / family$variance(g))
  aic = NULL
  if (!is.null(deviance) && !is.null(cumulative_weight)) {
    aic = family$aic(y, cumulative_weight, g, weights, deviance)    
  }
  RSS = sum(W*residuals^2)
  deviance = sum(family$dev.resids(y, g, weights))
  stop("here")
  null_dev = sum(family$dev.resids(y, , weights))
  list(XTWX=crossprod(mm, W * mm), XTWz=crossprod(mm, W*z), 
       deviance=deviance, null_dev=null_dev, cumulative_weight=nobs,
       aic=aic, RSS=RSS, contrasts=attr(mm, "contrasts"))
}

#' Data frame preprocessor 
#' 
#' @param col_types the column types of the data
#' @param col_names the column names of the data
#' @param sep the column separator
#' @param dfpp_fun the data.frame preprocessing function
#' @export
dfpp_gen = function(col_types, col_names=NULL, sep=",", 
                    dfpp_fun=function(x) x) {
  col_names = col_names
  col_types = col_types
  sep = sep
  dfpp_fun = dfpp_fun
  function(x) {
    x = dstrsplit(x, col_types=col_types, sep)
    colnames(x) = col_names
    dfpp_fun(x)
  }
}

#' Perform a generalized linear regression
#'
#' @param form an object of class ‘"formula"’ (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param data a connection to read data from
#' @param dfpp the data frame preprocessor. a function that turns the 
#' read data from a chunk into a properly formatted data frame
#' @param beta_start starting values for the linear predictor.
#' @param control a list of parameters for controlling the fitting process.
#' @param method the method to be used in fitting the model.  The default is
#' iteratively reweighted least squares (irls). Support for stochastic gradient
#' decent is coming.
#' @param contrasts an optional list. See the ‘contrasts.arg’ of ‘model.matrix.default’.
#' @param sep if using a connection or file, which character is used as a separator between elements?
#' @param parallel how many logical processor cores to use (default 1)
#' @export
ioglm = function(form, family = gaussian(), data, dfpp,
                 beta_start=NULL, control=list(maxit=25, epsilon=1e-08, 
                 trace=FALSE), method = "irls",
                 contrasts = NULL, sep=",", parallel=1, ...) {
  call = match.call()
  ret = NULL
  beta_old = beta_start
  cumulative_weight = NULL
  deviance = NULL
  if (is.data.frame(data)) {
    # The data.frame implementation.
    for (i in 1:control$maxit) {
      mm = model.matrix(form, data, contrasts)
      if (control$trace)
        cat("iteration", i, "\n")
      glm_m = glm_kernel(data[row.names(mm),all.vars(form)[1]], mm, family, 
                         beta_old, deviance, cumulative_weight)
      # TODO: Add checking for singularities here.
      XTWX = glm_m$XTWX
      XTWz = glm_m$XTWz
      deviance = glm_m$deviance
      aic = glm_m$aic
      RSS = glm_m$RSS
      null_dev=glm_m$null_dev
      cumulative_weight = glm_m$cumulative_weight
      contrasts = glm_m$contrasts
      beta = solve(XTWX, tol=2*.Machine$double.eps) %*% XTWz
      if (!is.null(beta_old) && 
          as.vector(sqrt(crossprod(beta-beta_old))) < control$epsilon) break 
      beta_old = beta
    }
  } else {
    # The iotools implementation.
    for (i in 1:control$maxit) {
      cvs = chunk.apply(data,
        function(x) {
          df = dfpp(x)
          mm = model.matrix(form, df, contrasts)
          glm_kernel(df[row.names(mm),all.vars(form)[1]], mm, family, beta_old,
                     deviance, cumulative_weight)
      }, CH.MERGE=list, parallel=parallel)
      XTWX = Reduce(`+`, Map(function(x) x$XTWX, cvs))
      XTWz = Reduce(`+`, Map(function(x) x$XTWz, cvs))
      aic = Reduce(`+`, Map(function(x) x$aic, cvs))
      RSS = Reduce(`+`, Map(function(x) x$RSS, cvs))
      null_dev = Reduce(`+`, Map(function(x) x$null_dev, cvs))
      contrasts = cvs[[1]]$contrasts
      deviance = Reduce(`+`, Map(function(x) x$deviance, cvs))
      cumulative_weight = Reduce(`+`, Map(function(x) x$cumulative_weight, cvs))
      # TODO: Add checking for singularities here.
      beta = solve(XTWX, tol=2*.Machine$double.eps) %*% XTWz
      if (!is.null(beta_old) && 
          as.vector(sqrt(crossprod(beta-beta_old))) < control$epsilon) break
      beta_old = beta
    }
  }
  if (as.vector(sqrt(crossprod(beta-beta_old))) < control$epsilon) {
    converged=TRUE
  }

  beta_names = row.names(beta)
  beta = as.vector(beta)
  names(beta) = beta_names

  rank=nrow(XTWX)

  aic_rest = ifelse(
    (family$family %in% c("Gamma", "inverse.gaussian", "gaussian")), 2, 0)
  aic = aic + 2 * nrow(XTWX) + aic_rest

  # This will have to change when we support weights.
  num_obs = cumulative_weight
  nulldf = num_obs - attributes(terms(form))$intercept
  resdf = num_obs - rank

  var_res = RSS/resdf
  dispersion = if (family$family %in% c("poisson", "binomial")) 1 else var_res

  if (missing(dfpp)) dfpp = NULL
  ret = list(coefficients=beta, 
    family=family,
    deviance = deviance,
    aic=aic,
    data=data,
    dfpp=dfpp,
    rank=rank,
    xtwx=XTWX,
    xtwz=XTWz,
    iter=i, 
    dispersion=dispersion,
    rss=RSS,
    converged=converged, 
    formula=form,
    call=call,
    num_obs=num_obs,
    nulldf=nulldf,
    null_dev=null_dev,
    resdf=resdf,
    terms=terms.formula(form),
    control=control,
    method=method,
    contrasts=contrasts)
  class(ret) = c("ioglm", "iolm")
  ret
}

#' Get the regression diagnostics for a linear regression
#' 
#' @param object an object return from ioglm
#' @param data a data.frame or connection to the data set where training was performed.
#' @param parallel how many logical processor cores to use (default 1)
#' @export
summary.ioglm = function(object, data, parallel=1, ...) {
  call = match.call()
  terms = object$terms
  if (missing(data)) {
    data = object$data
  }
  if (missing(parallel)) {
    parallel = object$parallel
  }
  dispersion = object$dispersion
  inv_scatter = solve(object$xtwx)
  standard_errors = sqrt(dispersion * diag(inv_scatter))
  stat_vals = object$coefficients/standard_errors
  if (object$family$family %in% c("binomial", "poisson")) {
    p_vals = 2 * pnorm(abs(stat_vals), lower.tail=FALSE)
  } else {
    p_vals = 2 * pt(abs(stat_vals), df=object$df, lower.tail=FALSE)
  }
  
  ret = list(call=call,
       terms=terms,
       family=object$family,
       deviance=object$deviance,
       aic=object$aic,
       contrasts=object$contrasts,
       df.residual=object$resdf,
 #      null.deviance=object$null_dev,
       df.null=object$nulldf,
       iter=object$iter,
       deviance.resid=NA,
       coefficients=object$coefficients,
       dispersion=object$dispersion,  
       df=c(ncol(object$xtwx), object$resdf, ncol(object$xtwx)),
       data=data,
       cov.unscaled=inv_scatter,
       cov.scaled=inv_scatter * dispersion)
  class(ret) = c("ioglm", "glm", "lm")
  ret
}
