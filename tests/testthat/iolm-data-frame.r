library(ioregression)
library(testthat)

form = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species

iofit = iolm(form, data=iris)
lmfit = lm(form, data=iris)

expect_equal(iofit$coefficients, lmfit$coefficients)
expect_equal(iofit$terms, lmfit$terms)
expect_equal(iofit$rank, lmfit$rank)
expect_equal(iofit$contrasts, lmfit$contrasts)

ios = summary(iofit)
lms = summary(lmfit)

expect_equal(ios$terms, lms$terms)
expect_equal(ios$coefficients, lms$coefficients, tol=0.001)
expect_equal(ios$aliased, lms$aliased)
expect_equal(ios$sigma, lms$sigma)
expect_equal(ios$df, lms$df)
expect_equal(ios$r.squared, lms$r.squared)
expect_equal(ios$adj.r.squared, lms$adj.r.squared)
expect_equal(ios$fstatistic, lms$fstatistic)
expect_equal(ios$cov.unscaled, lms$cov.unscaled)
