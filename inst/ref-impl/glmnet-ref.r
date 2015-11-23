
soft_thresh = function(x, g) {
  x = as.vector(x)
  w1 = which(g >= abs(x))
  w2 = which(g < abs(x) & x > 0)
  w3 = which(g < abs(x) & x < 0)
  ret = x
  ret[w1] = 0
  ret[w2] = x[w2]-g
  ret[w3] = x[w3]+g
  ret
}

lmnet_ref = function(X, y, lambda, alpha, maxit=10000, tol=1e-7) {
  beta = rep(1, ncol(X))
  xty = crossprod(X, y)
  xtx = crossprod(X)
  for(j in 1:maxit) {
    beta_old = beta
    beta = soft_thresh((xty - xtx %*% beta)/nrow(X) + beta, lambda*alpha)
    beta = beta / (1 + lambda*(1-alpha))
    if(sqrt(crossprod(beta-beta_old)) < tol) break
  }
  list(beta=beta, iterations=j)
}

x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit=glmnet(x,y, lambda=0.0454336178, standardize=FALSE, intercept=FALSE)

fit_ref = lmnet_ref(x, y, lambda=0.0454336178, alpha=1)
crossprod(fit_ref$beta, fit$beta)

max(abs(fit_ref$beta - fit$beta))

glmnet_ref = function(X, y, lambda, alpha, family=binomial, maxit=10, tol=1e-08)
{
  beta = matrix(rep(0,ncol(X)), ncol=1)
  resids = c()
  for(j in 1:maxit)
  {
    beta_outer_old = beta
    eta    = as.matrix(X %*% beta)
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (y - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    wx_norm = colSums(W*X^2)
    for (k in 1:maxit) {
      beta_inner_old = beta
      for (l in 1:length(beta)) {
      beta[l] = soft_thresh(sum(W*X[,l]*(z - X[,-l] %*% beta_inner_old[-l])), 
                              nrow(X)*lambda*alpha)
      }
      beta = beta / (wx_norm + lambda*(1-alpha))
      if(sqrt(as.double(crossprod(beta-beta_inner_old))) < tol) break
    }
    if (sqrt(as.double(crossprod(beta-beta_outer_old))) < tol) break
  }
  list(beta=beta,iterations=j)
}

x = matrix(rnorm(100*20),100,20) * -2
#x <- scale(x)
g2 = sample(0:1,100,replace=TRUE)

#lambda = 0
lambda = 0.0285513653

fg = glmnet(x, g2, family="gaussian", lambda=lambda, standardize=FALSE, intercept=FALSE)
fit = glmnet_ref(x, g2, family=gaussian, lambda=lambda, alpha=1, tol=0)

max(abs(fg$beta - fit$beta))
cbind(fg$beta, fit$beta)

fg = glmnet(x, g2, family="binomial", lambda=lambda, standardize=FALSE, intercept=FALSE)
fit = glmnet_ref(x, g2, family=binomial, lambda=lambda, alpha=1, tol=0)

print(max(abs(fg$beta - fit$beta)))
cbind(fg$beta, fit$beta)
stop("here")

fit$beta
fit$iterations
crossprod(fit$beta, fg$beta)

max(abs(fg$beta - fit$beta))



