
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

# Get xty, xtx_all, num_rows, and all_var_names from an adf.
net_matrices = function(data, formula, standardize, subset, weights, na.action,
                        offset, contrasts, parallel) {
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
    ret = list(xty=xty, xtx_all=xtx_all, num_rows=num_rows, 
               all_var_names=all_var_names)
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
      ret = list(xty=xty, xtx_all=xtx_all, num_rows=num_rows, 
                 all_var_names=all_var_names)
  }
  ret
}
