
#' Fit a generalized linear model with lasso or elasticnet regularization.
#'
#' Fit a generalized linear model via penalized maximum likelihood.
#' The regularization path is computed for the lasso or elasticnet
#' penalty at a grid of values for the regularization parameter
#' lambda. The function can deal data frames or abstract data frames.
#' @param formula the formula for the regression
#' @param family a description of the error distribution and link function
#' to be used in the model.
#' @param data an abstract data frame or something that can be coerced into one.
#' @param subset an options character string, which will be evaluated in the
#' frame of the data, to indicate which rows to include in the analysis
#' @param weights a optional character string, which will be evaluated in the
#' frame of the data, giving the sample weights for the
#' regression.
#' @param na.action a function which indicates what should happen when the data
#' contain 'NA's. See lm.fit for more details.
#' @param start starting values for the parameters in the linear predictor.
#' @param etastart starting values for the linear predictor.
#' @param mustart starting values for the vector of means.
#' @param offset a optional character string, which will be evaluated in the
#' frame of the data, giving the offsets for the regression
#' @param control a list of parameters for controlling the fitting process.
#' @param alpha the elasticnet mixing parameter 0 <= alpha <= 1.
#' @param lambda a user supplied value (or sequence of values) for the penalty
#' parameter. If not specified then a regularization path is created based
#' on the data.
#' @param contrasts contrasts to use with the regression. See the
#' \code{contrasts.arg} of \code{model.matrix.default}
#' @param standardize should the regression variables be normalized to have
#' mean zero and standard deviation one?
#' @param tol numeric tolerance.
#' @param max_it the maximum number of iterations per regression.
#' @param lambda_epsilon this value is multiplied by the max lambda in the
#' data to determine the minimum lambda when the regularization path is
#' determined from the data. This is ignored when lambda is specified.
#' @param nlambdas the number of lambdas to be generated in the regularization
#' path. This is ignored when lambda is specified.
#' @export
glmnet = function(formula, family, data, subset=NULL, weights=NULL, 
                   na.action=NULL, start=NULL, etastart, mustart,
                   offset=NULL, control=list(), alpha=1, 
                   lambda=NULL, contrasts=NULL, standardize=FALSE, tol=1e-7,
                   max_it=1e+05, lambda_epsilon=0.0001, nlambdas=100) {
  ret = ioirls(formula, family, data, weights, subset, na.action, start,
               etastart, mustart, offset, control, contrasts, trace,
               tol, glment_coordinate_descent_gen(lamda, alpha))
  class(ret) = c("ioglmnet")
  ret
}

partial_beta_kernel = function() {
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
  offset <- if (!is.null(d$offset)) d$offset else offset = rep.int(0,nobs)
  weights <- if (!is.null(d$w)) d$w else rep(1, nobs)
  family = passedVars$family

  eta = as.numeric(d$x %*% passedVars$beta)
  g = family$linkinv(eta <- eta + offset)
  gprime = family$mu.eta(eta)
  residuals = (d$y-g)/gprime
  z = eta - offset + residuals
  W = as.vector(weights * gprime^2 / family$variance(g))
  beta_new = rep(NA, length(beta))
  for (l in 1:length(beta_new)) {
    beta_new[l] = W * d$x[,l] * (z - d$x[,-l] %*% beta)
  }
  list(partial_beta=beta_new)
}

quad_loss_kernel = function() {
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
  offset <- if (!is.null(d$offset)) d$offset else offset = rep.int(0,nobs)
  weights <- if (!is.null(d$w)) d$w else rep(1, nobs)
  family = passedVars$family

  eta = as.numeric(d$x %*% passedVars$beta)
  g = family$linkinv(eta <- eta + offset)
  gprime = family$mu.eta(eta)
  residuals = (d$y-g)/gprime
  z = eta - offset + residuals
  W = as.vector(weights * gprime^2 / family$variance(g))
  list(partial_quad_loss = W * (z - d$x %*% beta)^2)
}

glmnet_coordinate_descent_gen(lambda, alpha) {
  function() {
    quad_loss_new = Inf
    for (i =1:control$maxit) {
      pvar = list(beta=beta, family=family)
      cvs = adf.apply(x=data, type="sparse.model", FUN=partial_beta_kernel,
                      args=pvar, formula=formula, subset=subset, 
                      weights=weights, na.action=na.action, offset=offset,
                      contrasts=contrasts)
      beta_new = Reduce(`+`, Map(function(x) x$partial_beta, cvs))        
      beta_new = soft_thresh(beta_new, cumulative_weight*lambda*alpha)
      beta_new = beta_new / (wx_norm + lambda *(1-alpha))

      pvar$beta = beta_new
      partial_quad_loss = adf.apply(x=data, type="sparse.model", 
                                    FUN=partial_beta_kernel,
                                    args=pvar, formula=formula, subset=subset,
                                    weights=weights, na.action=na.action, 
                                    offset=offset, contrasts=contrasts)
      quad_loss_new = Reduce(`+`, Map(function(x) x$partial_quad_loss, cvs))
      quad_loss_new = -1/2/nrow(X) * quad_loss_new +
        lambda * (1-alpha) * sum(beta_new^2)/2 + alpha * sum(beta_new)
      if (quad_loss > quad_loss_old) {
        beta = beta_new
        quad_loss = quad_loss_new
      }
      else break
    }
    beta
  }
}
