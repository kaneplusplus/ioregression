library(Matrix)
library(adf)
library(ioregression)
library(glmnet)

x=matrix(rnorm(100*20),100,20)
y=sample(1:2,100,replace=TRUE)
df = as.data.frame(x)
df$Y = y-1
names(df) 
names(df) = gsub("V", "X", names(df))

#alpha=0.5

fit=glmnet(x, y, family="binomial")
iofit = ioglmnet(Y ~ .-1, df, standardize=FALSE)


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
