library(Matrix)
library(glmnet)

data("QuickStartExample")

df = as.data.frame(x)
df$Y = as.vector(y)
names(df) = gsub("V", "X", names(df))

lambda=0.1
alpha=0.5

lmnet(Y ~ .-1, df, lambda=lambda, alpha=alpha)

col_means = colMeans(x)
x = x - matrix(data=col_means, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
col_sd = apply(x, 2, sd)
x = x / matrix(data=col_sd, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
y = (y-mean(y)) / sd(y)
glmnet(x, y, lambda = lambda, alpha=alpha)$beta
