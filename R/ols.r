
#' @param formula the formula for the regression
#' @param df the data frame to calculate xtx and xty for
#' @param contrasts the contrasts for categorical regressors
xtx_and_xty = function(formula, df, contrasts) {
  mm = model.matrix(formula, df, contrasts=contrasts)
  xtx = crossprod(mm)
  xty = crossprod(mm,
                  matrix(df[row.names(mm), as.character(formula)[2]],
                          ncol=1))
  sum_y = sum(df[row.names(mm),as.character(formula)[2]])
  n = nrow(mm)
  rm(df, mm)
  list(xtx=xtx, xty=xty, sum_y=sum_y, n=n)
}

# The chunk-wise operation for getting the RSS
yty_and_btxtxb= function(formula, b, df, contrasts, y_mean) {
  mm = model.matrix(formula, df, contrasts=contrasts)
  b_mat = matrix(b[colnames(mm)], nrow=ncol(mm), ncol=1)
  btxtxb = crossprod(mm %*% b_mat)
  yty = crossprod(matrix(df[row.names(mm), as.character(formula)[2]], ncol=1))
  syy = sum((df[row.names(mm), as.character(formula)[2]] - y_mean)^2)
  rss = sum((df[row.names(mm), as.character(formula)[2]] - mm %*% b_mat)^2)
  list(yty=yty, btxtxb=btxtxb, syy=syy, rss=rss)
}

#' Perform a big linear regression
#'
#' @param formula the formual for the regression
#' @param data a connection to read data from
#' @param data_frame_preprocessor a function that turns the data read from
#' a connection into a properly formatted data frame
#' @param contrasts the contrasts for categorical regressors
#' @param parallel the number of cores to apply to the regression
#' @export
blm = function(formula, data, data_frame_preprocessor=function(x) x, 
               contrasts=NULL, parallel=1, tol=-1) {
  call = match.call()
  if (is.data.frame(data)) {
    cvs = xtx_and_xty(formula, data_frame_preprocessor(data),contrasts)
    xtx = cvs$xtx
    xty = cvs$xty
    sum_y = cvs$sum_y
    n = cvs$n
  } else {
    if (!is.character(data)) {
      stop("Only data frames or file names are currently supported as input")
    }
    cvs = chunk.apply(iotools:::input(data),
      function(x) {
        xtx_and_xty(formula, data_frame_preprocessor(x), contrasts)
      },
      CH.MERGE=list, parallel=parallel)
    xtx = Reduce(`+`, Map(function(x) x$xtx, cvs))
    xty = Reduce(`+`, Map(function(x) x$xty, cvs))
    sum_y = Reduce(`+`, Map(function(x) x$sum_y, cvs))
    n = Reduce(`+`, Map(function(x) x$n, cvs))
  }
  design_matrix_names = colnames(xtx)
  # Get rid of colinear variables.
  ch = chol(xtx, pivot=TRUE, tol=tol)
  if (attr(ch, "rank") < ncol(xtx)) {
    # xtx is rank deficient. 
    new_name_order = attr(ch, "pivot")
    effective_rank = attr(ch, "rank")
    ft = terms(formula)
    formula_vars = colnames(attr(ft, "factors"))
    keep_vars = new_name_order[1:effective_rank]
    drop_var_names = colnames(xtx)[setdiff(new_name_order, keep_vars)]
    if (length(setdiff(drop_var_names, formula_vars)) > 0) {
      # A collinearity exists among factor levels in the design matrix.
      # Error out with an appropriate message.
      stop(paste("The following design matrix variables are colinear but",
                 "could not be dropped because they are associated with a",
                 "factor variable:",
                 setdiff(drop_var_names, formula_vars),
                 sep="\n"))
    } else {
      # The variable to drop can be dropped but we'll warn the user.
      warning(paste("Removing the following variables due to colinearity:",
              colnames(xtx)[tail(new_name_order,nrow(xtx)-effective_rank)],
              sep="\n"))
    } 
    #Reorder xtx and xty in terms of the pivots up to the rank.
    xtx = xtx[keep_vars, keep_vars]
    xty = xty[keep_vars,]
    # Now update the formula so that it doesn't have the colinear terms.
    drop_inds = match(drop_var_names, formula_vars)
    ft = drop.terms(ft, drop_inds, keep.response=TRUE)
    formula = formula(ft)
  }
  coefficients = as.numeric(solve(xtx) %*% xty)
  names(coefficients) = rownames(xtx)
  terms = terms(formula)
  contrasts = contrasts
  ret = list(coefficients=coefficients, call=call, terms=terms, 
             design_matrix_names=design_matrix_names, xtx=xtx, sum_y=sum_y,
             n=n)
  class(ret) = "blm"
  ret
}

#' Get the regression diagnostics for a big linear regression
#' 
#' @param object an object return from blm
#' @param data a data.frame or connection to the data set where training was performed.
#' @export
summary.blm = function(object, data, ...) {
  call = match.call()
  terms = object$terms
  if (is.data.frame(data)) {
    rss_chunk = yty_and_btxtxb(formula(object$terms), 
                               object$coefficients,
                               data_frame_preprocessor(data),
                               object$coefficients,
                               object$sum_y/object$n)
    yty = rss_chunk$yty
    btxtxb = rss_chunk$btxtxb
    syy = rss_chunk$syy
    rss = rss_chunk$rss
  } else {
    rss_chunks = chunk.apply(iotools:::input(data),
      function(x) {
        yty_and_btxtxb(formula(object$terms),
                       object$coefficients, 
                       data_frame_preprocessor(x), 
                       object$contrasts,
                       object$sum_y/object$n)
      },
      CH.MERGE=list, parallel=2)
    yty = Reduce(`+`, Map(function(x) x$yty, rss_chunks))
    btxtxb = Reduce(`+`, Map(function(x) x$btxtxb, rss_chunks))
    syy = Reduce(`+`, Map(function(x) x$syy, rss_chunks))
    rss = Reduce(`+`, Map(function(x) x$rss, rss_chunks))
  }
  rss_beta = yty - btxtxb
  df = c(length(object$design_matrix_names), 
         object$n-length(object$design_matrix_names),
         length(object$design_matrix_names))
  sigma = as.numeric(sqrt(rss_beta / df[2]))
  #vcv = sigma * sqrt(solve(object$xtx)) 
  cov_unscaled = solve(object$xtx)
  se = sigma * sqrt(diag(cov_unscaled))
  tv = object$coefficients/se
  ptv = apply(cbind( pt(tv, df), 1-pt(tv, df) ), 1, min)
  coefficients = cbind(object$coefficients, se, tv, ptv)
  colnames(coefficients) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  # Aliasing could be better.
  aliased = rep(FALSE, length(object$coefficients))
  names(aliased) = names(object$coefficients)
  r_squared = 1 - rss/syy
  adj_r_squared = 1 - (1-r_squared) * sum(df[1:2], -1)/df[2]
  fstatistic = c((syy - rss)/sigma^2/(df[1]-1), df[1]-1, dendf=df[2])
  names(fstatistic) = c("value", "numdf", "dendf")
  ret = list(call=call, terms=object$terms, coefficients=coefficients,
       aliased=aliased, sigma=sigma, df=df, r.squared=r_squared,
       adj.r.squared=adj_r_squared, fstatistic=fstatistic, 
       cov.unscaled=cov_unscaled)
  class(ret) = "summary.lm"
  ret
}
