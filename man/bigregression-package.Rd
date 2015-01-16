\name{bigregression-package}
\alias{bigregression-package}
\alias{bigregression}
\docType{package}
\title{
Out-of-core regressions with connection inputs
}
\description{
This package implements out-of-core regression routines (OLS, GLM,
LASSO, step-wise) and can be used with any input implemented as an R
connection.
}
\details{
\tabular{ll}{
Package: \tab bigregression\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-01-16\cr
License: \tab Apache License 2.0\cr
}

Index of functions:
\preformatted{
blm
}
}
\author{
Taylor Arnold, Michael J. Kane, Simon Urbanek

Maintainer: Michael J. Kane <kaneplusplus@gmail.com>
}
\examples{
blm(Sepal.Width ~ Sepal.Length + species, data=iris)
}
