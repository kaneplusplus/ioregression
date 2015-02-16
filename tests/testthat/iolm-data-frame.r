library(ioregression)
library(testthat)

form = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species
iofit = iolm(form, data=iris)

lmfit = lm(form, data=iris)

expect_equal(iofit$coefficients, lmfit$coefficients)
expect_equal(iofit$terms, lmfit$terms)
expect_equal(iofit$rank, lmfit$rank)
expect_equal(iofit$contrasts, lmfit$contrasts)
