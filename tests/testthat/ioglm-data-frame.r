library(ioregression)
library(testthat)

data = read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
form = admit ~ gre + gpa + rank
glmfit = glm(form, data=data, family = binomial())
iofit = ioglm(form, data = data, family = binomial())

expect_equal(iofit$coefficients, glmfit$coefficients)
expect_equal(iofit$rank, glmfit$rank)
expect_equal(iofit$family, glmfit$family)
expect_equal(iofit$deviance, glmfit$deviance)
expect_equal(iofit$aic, glmfit$aic)
expect_equal(iofit$iter, glmfit$iter+1)
expect_equal(iofit$converged, glmfit$converged)
expect_equal(iofit$formula, glmfit$formula)
expect_equal(iofit$terms, glmfit$terms)
expect_equal(iofit$contrasts, glmfit$contrasts)

sg = summary(glmfit)
s = summary(iofit)
