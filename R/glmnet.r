
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
#' @param lambda_epsilon this value is multiplied by the max lambda in the
#' data to determine the minimum lambda when the regularization path is
#' determined from the data. This is ignored when lambda is specified.
#' @param nlambdas the number of lambdas to be generated in the regularization
#' path. This is ignored when lambda is specified.
#' @export
glmnet = function(formula, family, data, subset=NULL, weights=NULL, 
                   na.action=NULL, offset=NULL, alpha=1, lambda=NULL, 
                   contrasts=NULL, standardize=FALSE, tol=1e-7,
                   max_it=1e+05, lambda_epsilon=0.0001, nlambdas=100) {

  # Under-development related messages.
  if (!missing(weights)) stop("Weights are not yet supported.")

  call = match.call()
  if (!inherits(data, "adf")) data = adf(data)

  if (!is.null(weights) && !is.character(weights <- weights[[1]]))
    stop("weights must be a length one character vector")
  if (!is.null(subset) && !is.character(subset <- subset[[1]]))
    stop("subset must be a length one character vector")
  if (!is.null(offset) && !is.character(offset <- offset[[1]]))
    stop("offset must be a length one character vector")
}
