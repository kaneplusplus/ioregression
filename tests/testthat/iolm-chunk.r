library(testthat)
library(ioregression)

# Download the data to a temp directory.
temp_dir = tempdir()
bz_file_name = file.path(temp_dir, "2008.csv.bz2")
file_name = file.path(temp_dir, "2008.csv")

if (!file.exists(bz_file_name)) {
    # Download the 2008 airline data.
    download.file("http://stat-computing.org/dataexpo/2009/2008.csv.bz2", 
                  bz_file_name)
}

# The data column names and types
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
form = ArrDelay ~ Distance 

# Use ioregression to create a linear model.
iofit = iolm(form, bzfile(bz_file_name, "rb"), 
             dfpp_gen(col_types, col_names, sep=","), parallel=4)

# Create a linear regression the R way.
x = read.csv(bzfile(bz_file_name, "r"), colClasses=col_types)
lmfit = lm(form, x)

# Check that the regression values are the same.
expect_equal(iofit$coefficients, lmfit$coefficients)
expect_equal(iofit$terms, lmfit$terms)
expect_equal(iofit$rank, lmfit$rank)
expect_equal(iofit$contrasts, lmfit$contrasts)

# Create summaries of the regression objects.
ios = summary(iofit, data=bzfile(bz_file_name, "rb"))
lms = summary(lmfit)

# Check that the summaries are the same.
expect_equal(ios$terms, lms$terms)
expect_equal(ios$coefficients, lms$coefficients, tol=0.001)
expect_equal(ios$aliased, lms$aliased)
expect_equal(ios$sigma, lms$sigma)
expect_equal(ios$df, lms$df)
expect_equal(ios$r.squared, lms$r.squared)
expect_equal(ios$adj.r.squared, lms$adj.r.squared)
expect_equal(ios$fstatistic, lms$fstatistic)
expect_equal(ios$cov.unscaled, lms$cov.unscaled)

