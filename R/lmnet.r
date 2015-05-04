
# The soft-threshold function
#
# @export
soft_thresh = function(x, gamma) {
  ret = sign(x) * (abs(x) - gamma)
  ret[abs(ret) < gamma] = 0
  ret
}

# The lmnet function
#
# @export
lmnet = function(formula, data, subset=NULL, weights=NULL, na.action=NULL,
                 offset=NULL, alpha=1, lambda=NULL, contrasts=NULL,
                 lambda.min.ratio=ifelse(nobs < nvars, 0.01, 1e-04),
                 standardize=TRUE, tolerance=1e-7, dfmax=nvars + 1, 
                 max_it = 1e+05, filter=c("strong", "safe", "none")) {

  # Under-development related messages.
  if (!standardize) stop("Unstandardizing data is not yet supported.")
  if (is.null(lambda)) stop("You still need to provide a value of lambda.")
  if (!missing(weights)) stop("Weights are not yet supported.")
  
  call = match.call()
  if (!inherits(data, "adf")) data = as.adf(data)

  # Get the standardization information as well as the
  # lambdas from the data.
  if (standardize) {
    stand_info = adf.apply(x=data, type="sparse.model",
      FUN=function(d,passedVars) {
        if (nrow(d$x) == 0L) return(NULL)
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
          sum_w = sum(d$w)
          d$x = d$x * sqrt(d$w)
          d$y = d$y * sqrt(d$w)
        } else {
          sum_w = nrow(d$x)
        }
        return(list(lambda = Matrix::crossprod(d$x, d$y),
                    n = nrow(d$x),
                    sum_w = sum_w,
                    sum_x = apply(d$x,2,sum),
                    sum_x_square = apply(d$x^2, 2, sum),
                    sum_y = sum(d$y * d$w)
                    sum_y_squared = sum((d$y * d$w)^2)
                    contrasts=attr(d$x, "contrasts")))

      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)
    stand_info = stand_info[!sapply(stand_info, is.null)]
    if (length(stand_info) == 0L) stop("No valid data.")
    data_lambdas = Reduce(`+`, Map(function(x) x$lambda, stand_info))    
    n = Reduce(`+`, Map(function(x) x$n , stand_info))    
    sum_w = Reduce(`+`, Map(function(x) x$sum_w , stand_info))    
    mean_x = Reduce(`+`, Map(function(x) x$sum_x, stand_info)) / n
    sum_x_square=Reduce(`+`, Map(function(x) x$sum_x_square, stand_info)) 
    sum_y = Reduce(`+`, Map(function(x) x$sum_y, stand_info)) 
    sum_y_squared = Reduce(`+`, Map(function(x) x$sum_y_squared, stand_info)) 
    contrasts=stand_info[[1]]$contrasts
  } else {
    stop("Non-standardized regressors are not yet supported.")
  }
  lambda_max = max(data_lambdas)

  # Rather than checking for a single lambda, this will be changed to 
  # iterate and create a regularization path.
  if (length(lambda) > 1L) stop("We don't support lambda trajectories yet.")
 
  # Next filter out columns
  if (filter[1] == "strong") {
    active_regressors = which(data_lambdas > 2*lambda - max(data_lambdas))
  } else if (filter[1] == "safe") {
    active_regressors = 
      which(data_lambdas > lambda - sqrt(sum_x_square)*sqrt(sum_y_square)*
        (lambda_max - lambda)/lambda_max)
  } else if (filter[1] == "none") {
    active_regressors = 1:length(data_lambdas)
  }
  else {
    stop("Unsupported filter type.")
  }

  # Now iterate until we get the slope coefficients. 
  # Note that for now I'm not going to remove unfiltered slope coefficients
  # that go to zero. This could be added.
  beta = Matrix(1, nrow=length(active_regressors))
  beta_old = -beta
  it_num = 0
  while (it_num <= max_it && as.vector(crossprod(beta-beta_old)) > tolerance) {
    beta_old = beta 

    # Pass through the data to get the covariance update information.
    cov_update_info = adf.apply(x=data, type="sparse.model",
      FUN=function(d,passedVars) {
        if (nrow(d$x) == 0L) return(NULL)
        if (!is.null(d$offset)) d$y = d$y - d$offset
        if (!is.null(d$w)) {
          if (any(d$w == 0)) {
            ok = d$w != 0
            d$w = d$w[ok]
            d$x = d$x[ok,,drop = FALSE]
            d$y = d$y[ok]
            if (!is.null(d$offset)) d$offset = d$offset[ok]
          }
          d$x = d$x * sqrt(d$w)
          d$y = d$y * sqrt(d$w)
        } else {
          sum_w = nrow(d$x)
        }
        if (standardize) {
          d$x = sqrt(
            (sum_x_square[active_regressors] - mean_x[active_regressors]^2)/n)*
            (d$x[,active_regressors] - mean_x[active_regressors])
        } else {
          # TODO: This should be easy to fix.
          stop("Unstandardized not supported")
        }
        return(list(xty = crossprod(d$x, d$y),
                    xtx = crossprod(d$x))
      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)
    cov_update_info = cov_update_info[!sapply(cov_update_info, is.null)]
    xty = Reduce(`+`, Map(function(x) x$xty, cov_update_info))
    xtx = Reduce(`+`, Map(function(x) x$xtx, cov_update_info))
    ud = xty - xtx %*% beta
    beta = soft_thresh(ud / n + beta, lambda*alpha) / (1 + lambda*(1-alpha))
    it_num = it_num + 1 
  }
  if (it_num > max_it)
    warning("The regression did not converge.")
  beta
}

