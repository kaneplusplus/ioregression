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
lambda = .1
alpha = 0.5
tol = 0.01

# elasticnet using no filtering of beta coefficients
while( as.vector(crossprod(beta-beta_old)) > tol ) {
  beta_old = beta
  beta = covariance_update(y, x, beta, lambda, alpha)
}

# Notice how close:
beta

# Is to:
glmnet(x, y, lambda = lambda, alpha=alpha)$beta


# Filter coefficients using the SAFE rules. Note this assumes 
# columns of x have been standardized.
# The functions returns whether or not you should keep beta at
# the corresponding positions (TRUE means keep it).
filter_safe = function(y, x, lambda, lambda_max=lambda) {
  xty = crossprod(x, y)
  xty >= lambda - sd(y) * (lambda_max - lambda) / lambda_max 
}

filter_strong = function(y, x, lambda, lambda_max=lambda) {
  crossprod(x, y) >= 2*lambda - lambda_max
}

# elasticnet using SAFE filtering of beta coefficients
beta = Matrix(as.numeric(filter_safe(y, x, lambda), nrow=ncol(x)))
beta_old = -beta
while( as.vector(crossprod(beta-beta_old)) > tol ) {
  beta_old = beta
  beta = covariance_update(y, x, beta, lambda, alpha)
}

# SAFE
beta

beta = Matrix(as.numeric(filter_strong(y, x, lambda), nrow=ncol(x)))
beta_old = -beta
while( as.vector(crossprod(beta-beta_old)) > tol ) {
  beta_old = beta
  beta = covariance_update(y, x, beta, lambda, alpha)
}

# STRONG
beta

# Is to:
glmnet(x, y, lambda = lambda, alpha=alpha)$beta

