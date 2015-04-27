library(Matrix)
library(glmnet)

data("QuickStartExample")

x=matrix(rnorm(100*20),100,20)
y=rnorm(100)

# Standardize.
col_means = colMeans(x)
x = x - matrix(data=col_means, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
col_sd = apply(x, 2, sd)
x = x / matrix(data=col_sd, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)

# The soft-threshold function
S = function(x, gamma) {
  ret = sign(x) * (abs(x) - gamma)
  ret[abs(ret) < gamma] = 0
  ret
}

# Covariance updating. Use if you can fit a pxp matrix in memory.
covariance_update = function(y, x, beta, lambda, alpha, tol=0.0001) {
  zero_inds = which(abs(beta) < tol )
  x[zero_inds,] = 0
  x[,zero_inds] = 0
  ud = crossprod(x, y) - crossprod(x) %*% beta
  S( ud/nrow(x) + beta, lambda*alpha ) / (1 + lambda*(1-alpha))
}

beta = Matrix(1, nrow=ncol(x))
beta_old = -beta
lambda = 0.5
alpha = 0.5
tol = 0.01

while( as.vector(crossprod(beta-beta_old)) > tol ) {
  beta_old = beta
  beta = covariance_update(y, x, beta, lambda, alpha)
}

# Notice how close:
beta

# Is to:
glmnet(x, y, lambda = seq(1, 0, by=-0.05), alpha=0.5)$beta[,"s10"]



