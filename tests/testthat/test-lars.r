library(testthat)
library(ioregression)
library(adf)

library(MASS)

# Download the data to a temp directory.
bz_file_name = file.path(temp_dir <- tempdir(), "1987.csv.bz2")
if (!file.exists(bz_file_name))
  download.file("http://stat-computing.org/dataexpo/2009/1987.csv.bz2", bz_file_name)

# Create an abstract data frame and a real data frame:
data = adf(bz_file_name, sep=",", header=TRUE, conMethod="bzfile")
data = allFactorLevels(data)
if (!exists('df')) df = read.table(bzfile(bz_file_name), header=TRUE, sep=",")

# Create an iofit and a local model matrix:
obj = iolm(DepDelay ~ Distance + Month + ArrTime - 1, data=data)
mf  =   lm(DepDelay ~ Distance + Month + ArrTime - 1, data=df, method="model.frame")
X   = model.matrix(DepDelay ~ Distance + Month + ArrTime - 1, mf)
y   = model.response(mf)

# Fit lars algorithms
iofit = iolm.lars(obj,  intercept=FALSE)
lmfit = lars::lars(X,y, intercept=FALSE)

# Compare methods:
expect_equal(iofit$lambda, lmfit$lambda)
expect_equal(iofit$df, lmfit$df)
expect_equal(iofit$R2, lmfit$R2)
expect_equal(iofit$RSS, lmfit$RSS)
expect_equal(iofit$Cp, lmfit$Cp)
expect_equal(as.numeric(iofit$beta), as.numeric(lmfit$beta))
