library(Matrix)
library(glmnet)

data("QuickStartExample")

df = as.data.frame(x)
df$Y = as.vector(y)
names(df) = gsub("V", "X", names(df))

lambda=0.1
alpha=0.5

lmnet(Y ~ .-1, df, lambda=lambda, alpha=alpha)
