
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

glmnet_ref = function(X, y, lambda, alpha, family=binomial, maxit=25, tol=1e-08)
{
  beta = rep(0,ncol(x))
  for(j in 1:maxit)
  {
    eta    = X %*% beta
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    WXy = W*X*(y - X %*% beta)
    WX2 = W*X^2
    beta_old = beta
    beta = soft_thresh(Matrix::colSums(WXy), lambda*alpha) /
           (Matrix::colSums(WX^2) + lambda*(1-alpha))
    if(sqrt(crossprod(beta-beta_old)) < tol) break
  }
  list(coefficients=beta,iterations=j)
}
