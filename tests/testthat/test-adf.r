library(testthat)
library(ioregression)

# Download the data to a temp directory.
bz_file_name = file.path(temp_dir <- tempdir(), "1987.csv.bz2")
if (!file.exists(bz_file_name))
  download.file("http://stat-computing.org/dataexpo/2009/1987.csv.bz2", bz_file_name)

# Create an abstract data frame and a real data frame:
data = adf(bz_file_name, sep=",", header=TRUE, conMethod="bzfile")
data = allFactorLevels(data)
df = read.table(bzfile(bz_file_name), header=TRUE, sep=",")

# Check data column names
expect_equal(data$colNames, names(df))

# Check column types
c0 = as.character(sapply(df,class))
c1 = as.character(data$colClasses)
c0[c0 == "factor"] = "character"
expect_equal(c0, c1)

# Check factor levels
expect_equal(data$levels$Origin, levels(df$Origin))
expect_equal(data$levels$Dest, levels(df$Dest))
expect_equal(data$levels$UniqueCarrier, levels(df$UniqueCarrier))