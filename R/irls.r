
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

  if (!is.null(weights) && !is.character(weights <- weights[[1]])) {
    stop("weights must be a length one character vector")
  }
  if (!is.null(subset) && !is.character(subset <- subset[[1]])) {
    stop("subset must be a length one character vector")
  }
  if (!is.null(offset) && !is.character(offset <- offset[[1]])) {
    stop("offset must be a length one character vector")
  }

  converged <- FALSE
  beta <- beta_old <- start
  wtdmu <- deviance <-  cumulative_weight <- NULL
  for (i in 1:control$maxit) {
    pvar <- list(beta = beta, cw = cumulative_weight, family = family,
                 deviance = deviance, wtdmu = wtdmu)
    cvs <- adf.apply(x = data, type="sparse.model", FUN = glm_kernel,
                     args = pvar, formula = formula, subset = subset,
                     weights = weights, na.action = na.action,
                     offset = offset, contrasts = contrasts)
    cvs <- cvs[!sapply(cvs, is.null)]
    num_rows <- Reduce(`+`, Map(function(x) x$num_rows, cvs))
    XTWX <- Reduce(`+`, Map(function(x) x$XTWX, cvs))
    XTWz <- Reduce(`+`, Map(function(x) x$XTWz, cvs))
    deviance <- Reduce(`+`, Map(function(x) x$deviance, cvs))
    cumulative_weight <- Reduce(`+`, Map(function(x) x$cumulative_weight, cvs))
    wtdmu <- if( "(Intercept)" %in% colnames(XTWX) ) {
      Reduce(`+`, Map(function(x) x$wy, cvs))/cumulative_weight
    } else {
      family$linkinv(0)
    }
    contrasts <- cvs[[1]]$contrasts
    wx_norm <- Reduce(`+`, Map(function(x) x$wx_norm, cvs))
    # TODO: Add checking for singularities here.
    beta <- beta_update_fun(XTWX, XTWz, tol)
    if (!is.null(beta_old)) {
      err <- as.vector(Matrix::crossprod(beta-beta_old))
    } else {
      err <- control$epsilon * 2
    }
    p_val <- Inf
    if ( (!is.null(beta_old) && (err < control$epsilon)) ) {
      converged <- TRUE
      break
    }
    if (trace) {
      cat(sprintf("Delta: %02.4f Deviance: %02.4f Iterations - %d\n",
                  err,deviance,i))
    }
    beta_old <- beta
  }

  # We only calculate these here, as we only care about the
  #  converging loop:
  aic <- Reduce(`+`, Map(function(x) x$aic, cvs))
  RSS <- Reduce(`+`, Map(function(x) x$RSS, cvs))
  nobs <- Reduce(`+`, Map(function(x) x$nobs, cvs))
  null_dev <- Reduce(`+`, Map(function(x) x$null_dev, cvs))
  rank <- nrow(XTWX)

  aic <- aic + 2 * nrow(XTWX)

  # This will have to change when we support weights.
  nulldf <- nobs - attributes(terms(formula))$intercept
  resdf <- nobs - rank

  var_res <- RSS/resdf
  dispersion <- if (family$family %in% c("poisson", "binomial")) 1 else var_res

  ret <- list(
    coefficients = beta,
    family = family,
    deviance = deviance,
    aic = aic,
    data = data,
    rank = rank,
    xtwx = XTWX,
    xtwz = XTWz,
    iter = i,
    dispersion = dispersion,
    rss = RSS,
    converged = converged,
    formula = formula,
    call = call,
    num_obs = nobs,
    df.null = nulldf,
    null.deviance = null_dev,
    df.residual = resdf,
    terms = terms.formula(formula),
    control = control,
    contrasts = contrasts)
  class(ret) <- c("ioglm", "iolm")
  ret
}

