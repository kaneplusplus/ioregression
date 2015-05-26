library(testthat)
library(ioregression)

library(MASS)

# Download the data to a temp directory.
bz_file_name = file.path(temp_dir <- tempdir(), "1987.csv.bz2")
if (!file.exists(bz_file_name))
  download.file("http://stat-computing.org/dataexpo/2009/1987.csv.bz2", bz_file_name)

# Create an abstract data frame and a real data frame:
data = adf(bz_file_name, sep=",", header=TRUE, conMethod="bzfile")
data = allFactorLevels(data)
df = read.table(bzfile(bz_file_name), header=TRUE, sep=",")

# Test on a simple regression
iolmObj = iolm(DepDelay ~ Distance + ArrDelay, data=data)
iofit = iolm.ridge(iolmObj, lambda=seq(0, 1, by = 0.1))
lmfit = lm.ridge(DepDelay ~ Distance + ArrDelay, df,
                  lambda=seq(0, 1, by = 0.1)*nrow(df))

iocoef = coef(iofit)
lmcoef = coef(lmfit)
attributes(lmcoef) = attributes(iocoef)
expect_equal(iocoef, lmcoef, tol=0.01)

# Test on a regression with factors
iolmObj = iolm(DepDelay ~ Distance + UniqueCarrier, data=data)
iofit = iolm.ridge(iolmObj, lambda=seq(0, 1, by = 0.1))
lmfit = lm.ridge(DepDelay ~ Distance + UniqueCarrier, df,
                  lambda=seq(0, 1, by = 0.1)*nrow(df))

iocoef = coef(iofit)
lmcoef = coef(lmfit)
attributes(lmcoef) = attributes(iocoef)
expect_equal(iocoef, lmcoef, tol=0.005)
