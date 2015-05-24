library(testthat)
library(ioregression)

TOL = 1e-4

# Download the data to a temp directory.
bz_file_name = file.path(temp_dir <- tempdir(), "1987.csv.bz2")
if (!file.exists(bz_file_name))
  download.file("http://stat-computing.org/dataexpo/2009/1987.csv.bz2",
                bz_file_name)

# Create an abstract data frame and a real data frame:
data = adf(bz_file_name, sep=",", header=TRUE, conMethod="bzfile")
data = allFactorLevels(data)
df = read.table(bzfile(bz_file_name), header=TRUE, sep=",")

# Run a basic glm model and check the terms
iofit = ioglm((DepDelay > 15) ~ Distance + UniqueCarrier, data=data,
              family=binomial)
lmfit =   glm((DepDelay > 15) ~ Distance + UniqueCarrier, data=df,
            family=binomial)
expect_equal(coef(iofit),coef(lmfit),tol=TOL)
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients,tol=TOL)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],tol=TOL)

# Run a basic glm model with a poisson loss function and check the terms
iofit = ioglm((DepDelay > 15) ~ Distance + UniqueCarrier, data=data,
              family=poisson)
lmfit =   glm((DepDelay > 15) ~ Distance + UniqueCarrier, data=df,
              family=poisson)
expect_equal(coef(iofit),coef(lmfit),tol=TOL)
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients,tol=TOL)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],tol=TOL)

# Run a glm model without an intercept (a lot of the
# metrics have special cases w/o an intercept) and check the terms
iofit = ioglm((DepDelay > 15) ~ Distance + UniqueCarrier - 1, data=data,
              family=binomial)
lmfit =   glm((DepDelay > 15) ~ Distance + UniqueCarrier - 1, data=df,
              family=binomial)
expect_equal(coef(iofit),coef(lmfit),tol=TOL)
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients, tol=TOL)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")], tol=TOL)

# Run a glm model, only using observations where DepDelay > 0
iofit = ioglm((DepDelay > 15) ~ Distance + UniqueCarrier, data=data,
              family=binomial, subset="DepDelay > 0")
lmfit =   glm((DepDelay > 15) ~ Distance + UniqueCarrier, data=df,
            family=binomial, subset=df$DepDelay > 0)
expect_equal(coef(iofit),coef(lmfit))
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")])

# Run a glm model with a weight term (use Month since
# it is non-negative and the best thing available).
iofit = ioglm((DepDelay > 15) ~ Distance + UniqueCarrier, data=data,
              family=binomial, weight="Month")
lmfit =   glm((DepDelay > 15) ~ Distance + UniqueCarrier, data=df,
              family=binomial, weight=df$Month)
expect_equal(coef(iofit),coef(lmfit),tol=TOL)
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients,tol=TOL)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")], tol=TOL)

# Run a glm model with an offset. Use ArrDelay.
iofit = ioglm((DepDelay > 15) ~ Distance + UniqueCarrier, data=data,
              family=binomial, offset="ArrDelay")
lmfit =   glm((DepDelay > 15) ~ Distance + UniqueCarrier, data=df,
              family=binomial, offset=df$ArrDelay)
expect_equal(coef(iofit),coef(lmfit),tol=TOL)
expect_equal(summary(iofit)$coefficients, summary(lmfit)$coefficients,tol=TOL)
expect_equal(summary(iofit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")],
             summary(lmfit)[c("sigma", "df", "r.squared", "adj.r.squared",
                              "fstatistic")], tol=TOL)
