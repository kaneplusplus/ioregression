
#' Fit a linear model with lasso or elasticnet regularization
#'
#' Fit a linear model via penalized maxiumum likelihood.
#' The regularization path is computed for the lasso or elasticnet
#' penalty at a grid of values for the regularization parameter
#' lambda. The function can deal data frames or abstract data frames.
#' @importFrom       Matrix crossprod colSums Matrix
#' @param formula the formula for the regression
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

  nm = net_matrices(data, formula, standardize, subset, weights, na.action,
                    offset, contrasts, parallel)
  # net_matrices returns xty, xtx, num_rows, and all_var_names. 
  attach(nm)

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
  # We are done with nm. Detach so we dont' get a warning the next time the
  # function is called.
  detach(nm)
  list(call=call, beta=beta_path, lambda = lambda_path)
}

