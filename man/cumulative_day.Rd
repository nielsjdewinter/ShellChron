% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumulative_day.r
\name{cumulative_day}
\alias{cumulative_day}
\title{Function to detect year transitions and calculate cumulative age of model results}
\usage{
cumulative_day(
  resultarray,
  plotyearmarkers = TRUE,
  export_peakid = TRUE,
  path = tempdir()
)
}
\arguments{
\item{resultarray}{Array containing the full results of
the optimized growth model}

\item{plotyearmarkers}{Should the location of identified year
transitions be plotted? \code{TRUE/FALSE}}

\item{export_peakid}{Should the result of peak identification
be plotted? \code{TRUE/FALSE}}

\item{path}{Export path (defaults to tempdir())}
}
\value{
A new version of the Julian Day tab of the resultarray
with Julian Day model estimates replaced by estimates of
cumulative age of the record in days.
}
\description{
Takes the result of iterative growth modelling and
transforms data from Julian Day (0 - 365) to cumulative
day of the shell age by detecting where transitions
from one year to the next occur and adding full years
(365 days) to simulations in later years.
}
\examples{
testarray <- array(NA, dim = c(40, 36, 9)) # Create empty array
# with correct third dimension
windowfill <- seq(50, 500, 50) \%\% 365 # Create dummy simulation data 
# (ages) to copy through the array
for(i in 6:length(testarray[1, , 1])){
    testarray[, i, 3] <- c(windowfill, rep(NA, length(testarray[, 1, 3]) - 
        length(windowfill)))
    windowfill <- c(NA, (windowfill + 51) \%\% 365)
}
# Add dummy /code{D} column.
testarray[, 1, 3] <- seq(1, length(testarray[, 1, 3]), 1)
# Add dummy YEARMARKER column
testarray[, 3, 3] <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
# Add dummy d18Oc column
testarray[, 2, 3] <- sin((2 * pi * (testarray[, 1, 3] - 8 + 7 / 4)) / 7)
testarray2 <- suppressWarnings(cumulative_day(testarray, FALSE, FALSE, tempdir()))
# Apply function on array
}
\references{
package dependencies: zoo 1.8.7; scales 1.1.0; graphics
function dependencies: peakid
}
