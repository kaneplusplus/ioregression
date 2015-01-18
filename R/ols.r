
#' @param formula the formula for the regression
#' @param df the data frame to calculate xtx and xty for
#' @param contrasts the contrasts for categorical regressors
xtx_and_xty = function(formula, df, contrasts) {
  mm = model.matrix(formula, df, contrasts=contrasts)
  xtx = crossprod(mm)
  xty = crossprod(mm,
                  matrix(df[row.names(mm), as.character(formula)[2]],
                          ncol=1))
  rm(df, mm)
  list(xtx=xtx, xty=xty)
}

# The chunk-wise operation for getting the RSS
yty_and_btxtxb = function(formula, b, df, contrasts) {
  mm = model.matrix(formula, df, contrasts=contrasts)
  b_mat = matrix(b[rownames(mm)], nrow=nrow(xtx), nrow=1)
  btxtxb = crossprod(mm %*% b_mat)
  yty = crossprod(matrix(df[row.names(mm), as.character(formula)[2]], ncol=1))
  list(yty=yty, btxtxb=btxtxb)
}

#' @param formula the formual for the regression
#' @param data_con a connection to read data from
#' @param data_frame_preprocessor a function that turns the data read from
#' a connection into a properly formatted data frame
#' @param contrasts the contrasts for categorical regressors
#' @param parallel the number of cores to apply to the regression
blm = function(formula, data_con, data_frame_preprocessor=function(x) x, 
               contrasts=NULL, parallel=1, tol=-1) {
  call = match.call()
  if (is.data.frame(data_con)) {
    xtx_and_xty(formula, data_frame_preprocessor(data_con), contrasts)
  } else {
    if (!is.character(data_con)) {
      stop("Only data frames or file names are currently supported as input")
    }
    cvs = chunk.apply(iotools:::input(data_con),
      function(x) {
        xtx_and_xty(formula, data_frame_preprocessor(x), contrasts)
      },
      CH.MERGE=list, parallel=parallel)
    xtx = Reduce(`+`, Map(function(x) x$xtx, cvs))
    xty = Reduce(`+`, Map(function(x) x$xty, cvs))
  }
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
                 setdiff(drop_var_names, formula_vars)
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
  ret = list(coefficients=coefficients, call=call, terms=terms)
  class(ret) = "blm"
  ret
}

# TODO: Add parallel option
summary.blm = function(object, data_con, ...) {
  call = match.call()
  terms = object$terms
  rss_chunks = chunk.apply(iotools:::input(data_con),
    function(x) {
      yty_and_btxtxb(formula(object$terms), object$coefficients, 
                     data_frame_preprocessor(x), object$contrasts, 2)
    },
    CH.MERGE=list, parallel=2)
  yty = Reduce(`+`, Map(function(x) x$yty, rss_chunks))
  btxtxb = Reduce(`+`, Map(function(x) x$btxtxb, rss_chunks))
  rss = yty - btxtxb
}
