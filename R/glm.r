source("ols.r")

glm_kernel = function(y, mm, family, beta) {
  eta = mm %*% beta
  g = family$linkinv(eta)
  gprime = family$mu.eta(eta)
  z = eta + (y - g)/gprime
  W = as.vector(gprime^2 / family$variance(g))
  list(XTWX=crossprod(mm, W * mm), XTWz=crossprod(mm, W*z))
}

#' Data frame preprocessor 
#' 
#' @param col_types the column types of the data
#' @param sep the column separator
#' @param dfpp_fun the data.frame preprocessing function
#' @export
dfpp_gen = function(col_names, col_types, sep, dfpp_fun=function(x) x) {
  col_names = col_names
  col_types = col_types
  sep = sep
  dfpp_fun = dfpp_fun
  function(x) {
    dfpp_fun(dstrsplit(x, col_types=col_types, sep))
  }
}


#' Perform a linear regression
#'
#' @param formula an object of class ‘"formula"’ (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param data a connection to read data from
#' @param dfpp the data frame preprocessor. a function that turns the 
#' read data from a chunk into a properly formatted data frame
#' @param beta_start starting values for the linear predictor.
#' @param maxit the maximum number of iterations.
#' @param method the method to be used in fitting the model.  The default is
#' iteratively reweighted least squares. Support for stochastic gradient
#' decent is coming.
#' @param contrasts an optional list. See the ‘contrasts.arg’ of ‘model.matrix.default’.
#' @param sep if using a connection or file, which character is used as a separator between elements?
#' @param parallel how many logical processor cores to use (default 1)
#' @param verbose show extra information.
#' @export
ioglm = function(formula, family = gaussian(), data, dfpp,
                 beta_start=NULL, maxit=25,
                 tol=1e-08, model = TRUE, method = "irls",
                 contrasts = NULL, sep=",", parallel=1, verbose=FALSE, ...) {
  ret = NULL
  beta_old = beta_start
  if (parallel != 1)
    warning("parallel argument is not working... setting parallel to 1.")
  if (is.data.frame(data)) {
    # The data.frame implementation.
    for (i in 1:maxit) {
      mm = model.matrix(form, data, contrasts)
      if (verbose)
        cat("iteration", i, "\n")
      if (is.null(beta_old)) beta_old = matrix(0, ncol=1, nrow=ncol(mm))
      glm_m = glm_kernel(data[row.names(mm),all.vars(form)[1]], mm, family, 
                         beta_old)
      # TODO: Add checking for singularities here.
      beta = solve(glm_m$XTWX, tol=2*.Machine$double.eps) %*% glm_m$XTWz
      if (as.vector(sqrt(crossprod(beta-beta_old))) < tol) break      
      beta_old = beta
    }
    ret = c(coefficients=beta, iterations=i)
  } else {
    if (!is.character(data))
      stop("Unsupported input type")
    if (missing(dfpp))
      stop("You must specify a data frame preprocessor function.")
    # The iotools implementation.
    for (i in 1:maxit) {
      cvs = chunk.apply(data,
        function(x) {
          df = dfpp(x)
          mm = model.matrix(form, df, contrasts)
          if (is.null(beta_old)) beta_old = matrix(0, ncol=1, nrow=ncol(mm))
          glm_kernel(df[row.names(mm),all.vars(form)[1]], mm, family, beta_old)
      }, CH.MERGE=list, parallel=parallel)
      XTWX = Reduce(`+`, Map(function(x) x$XTWX, cvs))
      XTWz = Reduce(`+`, Map(function(x) x$XTWz, cvs))
      # TODO: Add checking for singularities here.
      beta = solve(XTWX, tol=2*.Machine$double.eps) %*% XTWz
      if (!is.null(beta_old) && as.vector(sqrt(crossprod(beta-beta_old))) < tol)
        break
      beta_old = beta
    }
    ret = list(coefficients=beta, iterations=i)
  }
  class(ret) = "ioglm"
  ret
}

#' Get the regression diagnostics for a linear regression
#' 
#' @param object an object return from ioglm
#' @param data a data.frame or connection to the data set where training was performed.
#' @param data_frame_preprocessor any preprocessing that needs to be performed on the data
#' @param sep if using a connection or file, which character is used as a separator between elements?
#' @param parallel how many logical processor cores to use (default 1)
#' @export
summary.ioglm = function(object, data, data_frame_preprocessor=function(x) x,
                         sep=",", parallel=1, ...) {
}
