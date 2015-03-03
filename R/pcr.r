#' Recalculate iolm with Principal Component Regression
#'
#' This function takes a pre-computed iolm object,
#' and recalculates the regression vector with only
#' the top k-principal components.
#'
#' @param object     an iolm object with the desired formula
#' @param k          A scalar or vector of degrees to fit.
#' @param normalize  If TRUE, each variable is standardized to
#'                   have unit L2 norm, otherwise it is left alone.
#'                   Default is TRUE.
#' @export
iolm.pcr = function(object, k=seq(1,ncol(object$xtx)-attr(object$terms, "intercept")),
                    normalize=TRUE) {
  if (!inherits(object, "iolm"))
    stop("input to object must be an iolm object! See ?ioregression::iolm for more info.")

  n = object$n
  p = ncol(object$xtx)
  mean_y = object$sum_y / n
  mean_x = object$mean_x
  noms = names(out$coefficients)
  intercept = as.logical(attr(object$terms, "intercept"))

  k = unique(round(k))
  if (any(k <= 0L) || any(k > p-1))
    stop("Values of k must be between 1 and p-1.")

  XtY = as.matrix(object$xty)
  XtX = as.matrix(object$xtx)
  YtY = as.matrix(object$yty)

  # If the model matrix already has an intercept, remove it:
  if (intercept) {
    XtX = XtX[-1,-1]
    XtY = XtY[-1]
    mean_x = mean_x[-1]
    noms = noms[-1]
    p = p - 1
  } else mu = 0

  # Need to center X and Y for PCA:
  mu = mean_y
  XtY = XtY - n * mean_x * mean_y
  XtX = XtX - n * outer(mean_x,mean_x)

  if (normalize) {
    normx = sqrt(Matrix::diag(XtX))
    names(normx) <- NULL
    XtX = XtX / outer(normx,normx)
    XtY = XtY / normx
  } else normx = rep(1,p)

  ev = eigen(XtX)

  beta = matrix(0.0, ncol=p, nrow=length(k))

  for (j in 1:length(k)) {
    v = ev$vector[,1:k[j]]
    wtw = t(v) %*% XtX %*% v
    wty = t(v) %*% XtY
    beta[,j] = v %*% solve(wtw,wty)
  }

  beta = scale(t(beta), FALSE, normx)
  rownames(beta) = k
  MtM = XtX - outer(mean_x,mean_x)*n
  RSS = drop(YtY + Matrix::diag(beta %*% MtM %*% t(beta)) - 2 * beta %*% XtY
             - mu^2*n + 2 * n * mu * beta %*% mean_x)

  # Add offset back:
  beta = cbind(mu - beta %*% mean_x, beta)
  colnames(beta) = c("(Intercept)", noms)
  df = k

  out = list(coefficients=beta, df=df, RSS=RSS,
              iolm=out, call=match.call())
  class(out) = c("iolm.pcr")
  return(out)
}

#' Print an iolm.pcr object
#'
#' @method print iolm.pcr
#' @param x   an iolars object to print the summary of
#' @param ... other inputs; currently unused
#' @export
print.iolm.pcr = function (x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  mat = cbind(lambda=x$lambda, df=x$df, RSS=x$RSS, AIC=x$AIC,
              BIC=x$BIC)
  rownames(mat) = rep("",nrow(mat))
  print(mat)
}

