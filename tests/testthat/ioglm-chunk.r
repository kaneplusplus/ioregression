library(testthat)
library(ioregression)

temp_dir = tempdir()
bz_file_name = file.path(temp_dir, "2008.csv.bz2")
data_file_name = file.path(temp_dir, "2008.csv")

if (!file.exists(data_file_name)) {
    # Download the 2008 airline data.
    download.file("http://stat-computing.org/dataexpo/2009/2008.csv.bz2", 
                  bz_file_name)
    system(paste("bzip2 -d ", bz_file_name))
}

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

# Create the preprocessing function 
dfpp = dfpp_gen(col_names, col_types, sep=",", 
                function(x) { 
                  colnames(x) = col_names 
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


