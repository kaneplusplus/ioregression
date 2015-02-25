library(testthat)
library(ioregression)

# Download the data to a temp directory.
temp_dir = tempdir()
bz_file_name = file.path(temp_dir, "2008.csv.bz2")
data_file_name = file.path(temp_dir, "2008.csv")

if (!file.exists(data_file_name)) {
    # Download the 2008 airline data and decompress it.
    download.file("http://stat-computing.org/dataexpo/2009/2008.csv.bz2", 
                  bz_file_name)
    system(paste("bzip2 -d ", bz_file_name))
}

# The data column names and types.
col_names = c("Year", "Month", "DayofMonth", "DayOfWeek", "DepTime", 
              "CRSDepTime", "ArrTime", "CRSArrTime", "UniqueCarrier", 
              "FlightNum", "TailNum", "ActualElapsedTime", "CRSElapsedTime", 
              "AirTime", "ArrDelay", "DepDelay", "Origin", "Dest", "Distance", 
              "TaxiIn", "TaxiOut", "Cancelled", "CancellationCode", "Diverted", 
              "CarrierDelay", "WeatherDelay", "NASDelay", "SecurityDelay", 
              "LateAircraftDelay") 

col_types = c(rep("integer", 8), "character", "integer", "character", 
              rep("integer", 5), "character", "character", rep("integer", 4), 
              "character", rep("integer", 6)) 

# Here's the formula we'll fix on. 
form = Delayed ~ DayOfWeek 

# Create the preprocessing function.
dfpp = dfpp_gen(col_types, col_names, sep=",", 
                function(x) { 
                  # Get rid of the header row. 
                  if (x$UniqueCarrier[1] == "UniqueCarrier") 
                    x = x[-1,] 
                  x$Delayed = as.numeric(x$ArrDelay > 15) 
                  x$DayOfWeek = factor(x$DayOfWeek, 1:7, 
                                       labels=as.character(1:7)) 
                  x 
                })

# Compare this
iofit = ioglm(form, binomial(), data_file_name, dfpp)
print(iofit)

x = read.csv(data_file_name, as.is=TRUE)
x$Delayed = as.numeric(x$ArrDelay > 15)
x$DayOfWeek = factor(x$DayOfWeek, 1:7, labels=as.character(1:7))
glmfit = glm(form, family=binomial(), x)

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

expect_equal(s$terms, sg$terms)
expect_equal(s$family, sg$family)
expect_equal(s$deviance, sg$deviance)
expect_equal(s$aic, sg$aic)
expect_equal(s$contrasts, sg$contrasts)
expect_equal(s$df.residual, sg$df.residual)
expect_equal(s$null.deviance, sg$null.deviance)
expect_equal(s$df.null, sg$df.null)
expect_equal(s$coefficients, sg$coefficients)
expect_equal(s$dispersion, sg$dispersion)
expect_equal(s$df, sg$df)
expect_equal(s$cov.unscaled, sg$cov.unscaled, tol=0.001)
expect_equal(s$cov.scaled, sg$cov.scaled, tol=0.001)


