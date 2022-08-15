\name{SST_model}
\alias{SST Model}
\title{SST Model}
\usage{
SST_model()
}
\description{
Calculates the sea surface temperature based on d18O values from marine shells.
It is almost exactly the function `d18O_model()`, which calculates d18O values based on temperature values, mineralogy, and d18O values of the water, except this function does the opposite and calculates the temperature based on d18O values.
}
\examples{
# Create dummy d18O data
dist <- seq(1, 40, 1)
val <- sin((0.5 * pi * (dist)) / 5)+1

d18O <- cbind(dist, val)

# Run SST model function
SST <- SST_model(d18O, 0, "calcite")
}
