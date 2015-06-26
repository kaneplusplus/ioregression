
soft_thresh = function(x, g) {
  x = as.vector(x)
  w1 = which(g >= abs(x))
  w2 = which(g < abs(x) & x > 0)
  w3 = which(g < abs(x) & x < 0)
  ret = x
  ret[w1] = 0
  ret[w2] = x[w2]-g
  ret[w3] = x[w3]+g
  Matrix::Matrix(ret, nrow=length(x))
}

#' Fit a linear model with lasso or elasticnet regularization
#'
#' Fit a linear model via penalized maxiumum likelihood.
#' The regularization path is computed for the lasso or elasticnet
#' penalty at a grid of values for the regularization parameter
#' lambda. Can deal data frames or abstract data frames.
#' @importFrom       Matrix crossprod colSums Matrix
#' @param formula the formulat for the regression
#' @param data an abstract data frame or something that can be coerced into one.
#' @param subset an options character string, which will be evaluated in the
#' frame of the data, to indicate which rows to include in the analysis
#' @param weights a optional character string, which will be evaluated in the
#' frame of the data, giving the sample weights for the
#' regression.
#' @param na.action a function which indicates what should happen when the data
#' contain 'NA's. See lm.fit for more details.
#' @param offset a optional character string, which will be evaluated in the
#' frame of the data, giving the offsets for the regression
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
#' @param filter should filtering rules be used to remove variables that will
#' not be needed in the regression? Default is strong corresponding to
#' http://www-stat.stanford.edu/~tibs/ftp/strong.pdf.
#' @param lambda_epsilon this value is multiplied by the max lambda in the
#' data to determine the minimum lambda when the regularization path is
#' determined from the data. This is ignored when lambda is specified.
#' @param nlambdas the number of lambdas to be generated in the regularization
#' path. This is ignored when lambda is specified.
#' @param parallel    integer. the number of parallel processes to use in the
#'                     calculation (*nix only).
#' @export
iolmnet = function(formula, data, subset=NULL, weights=NULL, na.action=NULL,
                   offset=NULL, alpha=1, lambda=NULL, contrasts=NULL,
                   standardize=FALSE, tol=1e-7,
                   max_it = 1e+05, filter=c("strong", "safe", "none"),
                   lambda_epsilon=0.0001, nlambdas=100,parallel=1L) {

  # Under-development related messages.
  if (!missing(weights)) stop("Weights are not yet supported.")

  call = match.call()
  if (!inherits(data, "adf")) data = as.adf(data)

  if (!is.null(weights) && !is.character(weights <- weights[[1]]))
    stop("weights must be a length one character vector")
  if (!is.null(subset) && !is.character(subset <- subset[[1]]))
    stop("subset must be a length one character vector")
  if (!is.null(offset) && !is.character(offset <- offset[[1]]))
    stop("offset must be a length one character vector")

  # Here's where we get the initial xtx, xty and number of rows. Note that
  # if we standardize that we need an extra pass through the data; one
  # to get the mans, one to get the standardized inner-products. If you try
  # to do it in a single pass you will lose stability.
  if (standardize) {
    stand_info = adf.apply(x=data, type="sparse.model",
      FUN=function(d,passedVars) {
        if (nrow(d$x) == 0L) return(NULL)
        if (colnames(d$x)[1] == "(Intercept)") {
          # Get rid of the intercept column.
          d$x = d$x[,2:ncol(d$x),drop=FALSE]
        }
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
        return(list(num_rows = nrow(d$x),
                    sum_w = sum_w,
                    col_sum_x = Matrix::colSums(d$x),
                    col_sum_x_squared = Matrix::colSums(d$x^2),
                    sum_y = sum(d$y),
                    contrasts=attr(d$x, "contrasts"),
                    all_var_names = colnames(d$x), 
                    warn_intercept=attributes(d$mt)$intercept))

      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts,
        parallel=parallel)
    stand_info = stand_info[!sapply(stand_info, is.null)]
    if (length(stand_info) == 0L) stop("No valid data.")
    num_rows = Reduce(`+`, Map(function(x) x$num_rows, stand_info))    
    sum_w = Reduce(`+`, Map(function(x) x$sum_w , stand_info))    
    col_sum_x = Reduce(`+`, Map(function(x) x$col_sum_x, stand_info))
    col_mean_x = col_sum_x / num_rows
    sum_y = Reduce(`+`, Map(function(x) x$sum_y, stand_info)) 
    mean_y = sum_y / num_rows
    col_sum_x_squared = Reduce(`+`,
      Map(function(x) x$col_sum_x_squared, stand_info))
    contrasts=stand_info[[1]]$contrasts
    all_var_names = stand_info[[1]]$all_var_names
    if (stand_info[[1]]$warn_intercept) 
      warning("Intercept is ignored when variables are standardized")

    # Get the standard deviations and then fix the lambdas.
    stand_info = adf.apply(x=data, type="sparse.model",
      FUN=function(d,passedVars) {
        if (nrow(d$x) == 0L) return(NULL)
        if (colnames(d$x)[1] == "(Intercept)") {
          # Get rid of the intercept column.
          d$x = d$x[,2:ncol(d$x),drop=FALSE]
        }
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
        x_bar = Matrix::Matrix(col_mean_x,ncol=ncol(d$x), nrow=nrow(d$x), 
                               byrow=TRUE)
        x_centered=d$x-x_bar
        y_centered=d$y-mean_y
        return(list(col_x_square_diff=Matrix::colSums(x_centered^2),
                    y_square_diff=sum(y_centered^2),
                    centered_xty=Matrix::crossprod(x_centered, y_centered),
                    all_var_names=colnames(d$x)))

      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)
      col_x_square_diff =
        Reduce(`+`, Map(function(x) x$col_x_square_diff, stand_info))
      x_sd = sqrt(col_x_square_diff / (num_rows-1))
      y_square_diff = Reduce(`+`, Map(function(x) x$y_square_diff, stand_info))
      y_sd = sqrt(y_square_diff / (num_rows-1))
      xty = Reduce(`+`,
        Map(function(x) x$centered_xty, stand_info)) / x_sd / y_sd
      all_var_names = stand_info[[1]]$all_var_names

    # Second pass in the standardized case.
    cov_update_info = adf.apply(x=data, type="sparse.model",
      FUN=function(d,passedVars) {
        #d$x = d$x[,active_regressors, drop=FALSE]
        if (nrow(d$x) == 0L) return(NULL)
        if (colnames(d$x)[1] == "(Intercept)") {
          # Get rid of the intercept column.
          d$x = d$x[,2:ncol(d$x),drop=FALSE]
        }
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
        d$x=(d$x-Matrix::Matrix(col_mean_x,
             ncol=ncol(d$x), nrow=nrow(d$x), byrow=TRUE)) /
             Matrix::Matrix(x_sd, ncol=ncol(d$x), nrow=nrow(d$x),
                    byrow=TRUE)
        d$y = (d$y - mean_y) / y_sd
        return(list(xtx = Matrix::crossprod(d$x)))
      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)
    cov_update_info = cov_update_info[!sapply(cov_update_info, is.null)]
    xtx_all = Reduce(`+`, Map(function(x) x$xtx, cov_update_info))
  } else {

    # Unstandardized.
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
          d$x = d$x * sqrt(d$w)
          d$y = d$y * sqrt(d$w)
        } else {
          sum_w = nrow(d$x)
        }
        return(list(num_rows=nrow(d$x), xty=Matrix::crossprod(d$x, d$y),
                    xtx=Matrix::crossprod(d$x),
                    contrasts=attr(d$x, "contrasts"),
                    all_var_names=colnames(d$x)))

      },formula=formula,subset=subset,weights=weights,
        na.action=na.action, offset=offset, contrasts=contrasts)
      xty = Reduce(`+`, Map(function(x) x$xty, stand_info))
      xtx_all = Reduce(`+`, Map(function(x) x$xtx, stand_info))
      num_rows = Reduce(`+`, Map(function(x) x$num_rows, stand_info))    
      all_var_names=stand_info[[1]]$all_var_names
  }

  data_lambdas = abs(xty) / num_rows / alpha
  lambda_k = max(abs(data_lambdas))
  if (is.null(lambda)) {
    lambda_path = seq(from=lambda_k, to=lambda_epsilon*lambda_k,
                      length.out=nlambdas)
  } else {
    lambda_path = lambda
  }
  beta_path = NULL
  lambda = c()

  # Note that we'll need to do something with the a0 when we're not
  # standardizing. For now they are 0.
  # The following could be parallelized.
  for(lambda in lambda_path) {
    if (filter[1] == "strong") {
      active_regressors =
        which(as.vector(data_lambdas) > 2*lambda - lambda_k)
    } else if (filter[1] == "safe") {
      active_regressors =
        which(as.vector(data_lambdas) >
          lambda - sqrt(sum(col_x_square_diff))*sqrt(sum(y_square_diff))*
          (lambda_k - lambda)/lambda_k)
    } else if (filter[1] == "none") {
      active_regressors = 1:length(data_lambdas)
    } else {
      stop("Unsupported filter type.")
    }
    xtx = xtx_all[active_regressors,,drop=FALSE]
    xtx = xtx[,active_regressors,drop=FALSE]
    if (nrow(xtx) > 0) {
      beta = Matrix::Matrix(1, nrow=nrow(xtx), ncol=1)
      beta_old = -beta
      it_num = 0
      while (it_num <= max_it &&
             as.vector(Matrix::crossprod(beta-beta_old)) > tol) {
        beta_old = beta
        ud = xty[active_regressors,,drop=FALSE] - xtx %*% beta
        beta = soft_thresh(ud / num_rows + beta, lambda*alpha) /
          (1 + lambda*(1-alpha))
        it_num = it_num + 1
        if (sum(beta == 0) >= length(beta) / 2)
          beta = suppressWarnings(as(beta, "dgCMatrix"))
      }
      if (it_num > max_it)
        warning("The regression did not converge.")
      beta_ret = Matrix::Matrix(0, ncol=1, nrow=(length(all_var_names)),
        dimnames=list(all_var_names, NULL))
      beta_ret[colnames(xtx),1] = beta
    } else {
      beta_ret = Matrix::Matrix(0, ncol=1, nrow=(length(all_var_names)),
        dimnames=list(all_var_names, NULL))
    }
    if (is.null(beta_path)) {
      beta_path = beta_ret
    } else {
      beta_path = cbind(beta_path, beta_ret)
    }
  }
  # TODO: this needs to be fixed.
  a0 = rep(0, length(lambda_path)) 
  list(call=call, a0=a0, beta=beta_path, lambda = lambda_path)
}

