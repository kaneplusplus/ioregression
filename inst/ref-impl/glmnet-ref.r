
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
  beta = rep(0, ncol(X))
  xty = crossprod(X, y)
  xtx = crossprod(X)
  for(j in 1:maxit) {
    beta_old = beta
    #beta = soft_thresh((xty - xtx %*% beta)/nrow(X) + beta, lambda*alpha)
    for (l in 1:length(beta)) {
      beta[l] = soft_thresh(sum(X[,l]*(y - X[,-l] %*% beta_old[-l])),
                            nrow(X)*lambda*alpha)
    }
    beta = beta / (colSums(X^2) + lambda*(1-alpha))
    if(sqrt(crossprod(beta-beta_old)) < tol) break
  }
  list(beta=beta, iterations=j)
}

set.seed(3)
x=matrix(rnorm(100*20),100,20) * 100
y= sample(0:1,100,replace=TRUE)
fit=glmnet(x,y, lambda=0.0454336178, standardize=FALSE, intercept=FALSE)

fit_ref = lmnet_ref(x, y, lambda=0.0454336178, alpha=1)
crossprod(fit_ref$beta, fit$beta)
print("LMNET")

print(cbind(fit$beta, fit_ref$beta))
stop("here")
max(abs(fit_ref$beta - fit$beta))

glmnet_ref = function(X, y, lambda, alpha, family=binomial, maxit=10, tol=1e-08)
{
  beta = matrix(rep(0,ncol(X)), ncol=1)
  ql = c()
  for(j in 1:maxit)
  {
    beta_outer_old = beta
    eta    = as.matrix(X %*% beta)
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (y - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    wx_norm = colSums(W*X^2)
    quad_loss = Inf
    for (k in 1:maxit) {
      quad_loss_old = Inf
      beta_inner_old = beta
      for (l in 1:length(beta)) {
        beta[l] = soft_thresh(sum(W*X[,l]*(z - X[,-l] %*% beta_inner_old[-l])), 
                              sum(W)*lambda*alpha)
      }
      beta = beta / (wx_norm + lambda*(1-alpha))
      quad_loss = -1/2/nrow(X) * sum(W*(z - X %*% beta)^2) + 
        lambda * (1-alpha) * sum(beta^2)/2 + alpha * sum(beta)
      if (quad_loss > quad_loss_old) quad_loss_old = quad_loss
      else break
      #if(sqrt(as.double(crossprod(beta-beta_inner_old))) < tol) break
    }
    ql = c(ql, quad_loss)
    if (sqrt(as.double(crossprod(beta-beta_outer_old))) < tol) break
  }
  list(beta=beta,iterations=j, ql=ql)
}

set.seed(3)
x = matrix(rnorm(100*20),100,20) * 100 
#x <- scale(x)
g2 = sample(0:1,100,replace=TRUE)

#lambda = 0
lambda = 0.0454336178

fg = glmnet(x, g2, family="gaussian", lambda=lambda, standardize=FALSE, intercept=FALSE)
fit = glmnet_ref(x, g2, family=gaussian, lambda=lambda, alpha=1, tol=0)

max(abs(fg$beta - fit$beta))
print("GLMNET GAUSSIAN")
print(cbind(fg$beta, fit$beta))

fg = glmnet(x, g2, family="binomial", lambda=lambda, standardize=FALSE, intercept=FALSE)
fit = glmnet_ref(x, g2, family=binomial, lambda=lambda, alpha=1, tol=0)

print(max(abs(fg$beta - fit$beta)))
print("GLMNET BINOMIAL")
print(cbind(fg$beta, fit$beta))

fit$beta
fit$iterations
crossprod(fit$beta, fg$beta)

max(abs(fg$beta - fit$beta))



