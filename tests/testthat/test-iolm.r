library(testthat)
library(adf)
library(ioregression)

# Download the data to a temp directory.
library(adf)
bz_file_name = file.path(temp_dir <- tempdir(), "1987.csv.bz2")
if (!file.exists(bz_file_name)) {
  download.file("http://stat-computing.org/dataexpo/2009/1987.csv.bz2", 
                bz_file_name)
}

# Create an abstract data frame and a real data frame:
data = adf(bz_file_name, conMethod="bzfile", header=TRUE, sep=",")
data = allFactorLevels(data)
if(!exists("df")) df = read.table(bzfile(bz_file_name), header=TRUE, sep=",")

# Run a basic linear regression and check the terms
iofit = iolm(DepDelay ~ Distance + UniqueCarrier, data=data)
lmfit =   lm(DepDelay ~ Distance + UniqueCarrier, data=df)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma","df","r.squared","adj.r.squared","fstatistic")],
             summary(lmfit)[c("sigma","df","r.squared","adj.r.squared","fstatistic")])

# Run a linear regression without an intercept (a lot of the
# metrics have special cases w/o an intercept) and check the terms
iofit = iolm(DepDelay ~ Distance + UniqueCarrier - 1, data=data)
lmfit =   lm(DepDelay ~ Distance + UniqueCarrier - 1, data=df)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma","df","r.squared","adj.r.squared","fstatistic")],
             summary(lmfit)[c("sigma","df","r.squared","adj.r.squared","fstatistic")])

# Run a basic linear regression, only using observations where DepDelay > 0
iofit = iolm(DepDelay ~ Distance + UniqueCarrier, data=data, 
             subset="DepDelay > 0")
lmfit =   lm(DepDelay ~ Distance + UniqueCarrier, data=df,   
             subset=df$DepDelay > 0)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")])

# Run a linear regression with a weight term (use Month since
# it is non-negative and the best thing available).
iofit = iolm(DepDelay ~ Distance + UniqueCarrier, data=data, weight="Month")
lmfit =   lm(DepDelay ~ Distance + UniqueCarrier, data=df,  weight=df$Month)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")])

# Run a linear regression with an offset. Use ArrDelay.
iofit = iolm(DepDelay ~ Distance + UniqueCarrier, data=data, offset="ArrDelay")
lmfit =   lm(DepDelay ~ Distance + UniqueCarrier, data=df,  offset=df$ArrDelay)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")])

# Test data frames
iofit = iolm(Sepal.Width ~ Sepal.Length, iris)
lmfit =   lm(Sepal.Width ~ Sepal.Length, iris)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")])

# Test dot notation
iofit = iolm(Sepal.Width ~ ., iris)
lmfit =   lm(Sepal.Width ~ ., iris)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")])
