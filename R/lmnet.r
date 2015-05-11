
# Soft-threshold a vector of values
#
# @export
soft_thresh = function(x, g) {
  x = as.vector(x)
  w1 = which(g >= abs(x))
  w2 = which(g < abs(x) & x > 0)
  w3 = which(g < abs(x) & x < 0)
  ret = x
  ret[w1] = 0
  ret[w2] = x[w2]-g
  ret[w3] = x[w3]+g
  Matrix(ret, nrow=length(x))
}

# Fit the lmnet 
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
  # lambdas from the data. Note that this will take two passes. One to get
  # the means, one to get the variances. If you try to do it in a single
  # pass you start to lose stability.
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
        return(list(n = nrow(d$x),
                    sum_w = sum_w,
                    sum_x = Matrix::colSums(d$x),
                    sum_x_squared = Matrix::colSums(d$x^2),
                    sum_y = sum(d$y * d$w),
                    sum_y_squared = sum((d$y * d$w)^2),
                    contrasts=attr(d$x, "contrasts"),
                    all_var_names = colnames(d$x)))

      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)

    stand_info = stand_info[!sapply(stand_info, is.null)]
    if (length(stand_info) == 0L) stop("No valid data.")
    n = Reduce(`+`, Map(function(x) x$n , stand_info))    
    sum_w = Reduce(`+`, Map(function(x) x$sum_w , stand_info))    
    mean_x = Reduce(`+`, Map(function(x) x$sum_x, stand_info)) / n
    sum_y = Reduce(`+`, Map(function(x) x$sum_y, stand_info)) 
    sum_y_squared = Reduce(`+`, Map(function(x) x$sum_y_squared, stand_info)) 
    sum_x_squared = Reduce(`+`, Map(function(x) x$sum_x_squared, stand_info)) 
    contrasts=stand_info[[1]]$contrasts
    all_var_names = stand_info[[1]]$all_var_names

    # Get the standard deviations and then fix the lambdas.
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
        x_centered=d$x-Matrix(mean_x, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
        return(list(x=d$x,
          square_diff= Matrix::colSums(x_centered^2),
          data_lambda_unnormalized=Matrix::crossprod(x_centered, d$y)))

      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)
      square_diff = Reduce(`+`, Map(function(x) x$square_diff, stand_info))
      x_sd = sqrt(square_diff / (n-1))
      data_lambdas = Reduce(`+`, Map(function(x) x$data_lambda_unnormalized,
                       stand_info)) / x_sd
      xty = data_lambdas

  } else {
    stop("Non-standardized regressors are not yet supported.")
  }
  lambda_max = max(data_lambdas)

  # Rather than checking for a single lambda, this will be changed to 
  # iterate and create a regularization path.
  if (length(lambda) > 1L) stop("We don't support regularization paths yet.")
 
  # Next filter out columns. Note that the active_regressor variable will
  # keep track of which variables go into the covariance updating.
  if (filter[1] == "strong") {
    active_regressors = 
      which(as.vector(data_lambdas) > 2*lambda - max(data_lambdas))
  } else if (filter[1] == "safe") {
    # TODO: we need ||x||_2 and ||y||_2 for this.
    active_regressors = 
      which(as.vector(data_lambdas) > 
        lambda - sqrt(sum(sum_x_squared))*sqrt(sum(sum_y_squared))*
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
  # Pass through the data to get the covariance update information.
  cov_update_info = adf.apply(x=data, type="sparse.model",
    FUN=function(d,passedVars) {
      d$x = d$x[,active_regressors]
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
        # Center
        # TODO: pick out the active regressors. Zero the rest.
        d$x=(d$x-Matrix(mean_x, ncol=ncol(d$x), nrow=nrow(d$x), byrow=TRUE)) /
          Matrix(x_sd, ncol=ncol(d$x), nrow=nrow(d$x), byrow=TRUE)
      } else {
        # TODO: This should be easy to fix.
        stop("Unstandardized not supported")
      }
      return(list(xtx = Matrix::crossprod(d$x)))
    },formula=formula,subset=subset,weights=weights,
      na.action=na.action, offset=offset, contrasts=contrasts)
  cov_update_info = cov_update_info[!sapply(cov_update_info, is.null)]
  
  xtx = Reduce(`+`, Map(function(x) x$xtx, cov_update_info))
  beta = Matrix(1, nrow=nrow(xtx))
  beta_old = -beta
  it_num = 0
  while (it_num <= max_it && 
         as.vector(Matrix::crossprod(beta-beta_old)) > tolerance) {
    beta_old = beta 
    ud = xty - xtx %*% beta
    beta = soft_thresh(ud / n + beta, lambda*alpha) / (1 + lambda*(1-alpha))
    it_num = it_num + 1 
    if (sum(beta == 0) >= length(beta) / 2) 
      beta = suppressWarnings(as(beta, "dgCMatrix"))
  }
  if (it_num > max_it)
    warning("The regression did not converge.")
  list(beta=beta, cn=colnames(xtx))
  beta_ret = Matrix(0, ncol=1, nrow=(length(all_var_names)), 
    dimnames=list(all_var_names, NULL))
  beta_ret[colnames(xtx),1] = beta
  beta_ret
}

