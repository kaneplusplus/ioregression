library(Matrix)
library(glmnet)

data("QuickStartExample")

df = as.data.frame(x)
df$Y = as.vector(y)
names(df) = gsub("V", "X", names(df))

alpha=0.5

iofit = iolmnet(Y ~ .-1, df, alpha=alpha, standardize=FALSE)
fit = glmnet(x, y, alpha=alpha)

print(fit$beta[,ncol(fit$beta)])
print(iofit$beta[,ncol(iofit$beta)])

iofit = iolmnet(Y ~ ., df, alpha=alpha, standardize=TRUE)
col_means = colMeans(x)
x = x - matrix(data=col_means, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
col_sd = apply(x, 2, sd)
x = x / matrix(data=col_sd, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
y = (y-mean(y)) / sd(y)
fit = glmnet(x, y, alpha=alpha)

print(fit$beta[,ncol(fit$beta)])
print(iofit$beta[,ncol(iofit$beta)])
