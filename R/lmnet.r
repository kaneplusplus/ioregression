
# The soft-threshold function
#
# @export
softmax = function(x, gamma) {
  ret = sign(x) * (abs(x) - gamma)
  ret[abs(ret) < gamma] = 0
  ret
}

# The lmnet function
#
# @export
lmnet = function(formula, data, weights, offset=NULL, alpha=1,
  lambda=NULL, lambda.min.ratio=ifelse(nobs < nvars, 0.01, 1e-04),
  standardize=TRUE, thresh=1e-7, dfmax=nvars + 1, 
  maxit = 1e+05, filter=c("strong", "safe")) {

  # Under-development related messages.
  if (!standardize) stop("Unstandardizing data is not yet supported.")
  if (is.null(lambda)) stop("You still need to provide a value of lambda.")
  if (!missing(weights)) stop("Weights are not yet supported.")
  
  call = match.call()
  if (!inherits(data, "adf")) data = as.adf(data)

  
  if (filter == "strong") {
    cvs = adf.apply(x=data, type="sparse.model",
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
                sum_y = sum(d$y)
                sum_w = nrow(d$x)
              }
              return(list(xtx = Matrix::crossprod(d$x),
                          xty = Matrix::crossprod(d$x, d$y),
                          yty = Matrix::crossprod(d$y),
                          n = nrow(d$x),
                          sum_y = sum_y,
                          sum_w = sum_w,
                          mean_x = apply(d$x,2,sum),
                          contrasts=attr(d$x, "contrasts")))

            },formula=formula,subset=subset,weights=weights,
              na.action=na.action, offset=offset, contrasts=contrasts)
  }
  else {
    stop("Unsupported filter type.")
  }
}

