% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import_object.r
\name{data_import_object}
\alias{data_import_object}
\title{Function to prepare d18O data from a preexisting object
and process yearmarkers and calculation windows}
\usage{
data_import_object(input_object)
}
\arguments{
\item{input_object}{Name of the object that contains sampling distance and
d18O data. Note that sampling distance should be given in micrometers,
because the SCEUA model underperforms when the growth rate figures are very
small (<0.1 mm/day).}
}
\value{
A list containing an object with the original data and details on
the position and length of modeling windows
}
\description{
Takes the name of an object, reads the user-provided yearmarkers in the file
and uses them as a basis for the length of windows used throughout the model.
This ensures that windows are not too short and by default contain at least
one year of growth for modeling.
}
\examples{
Create virtual data
dat <- as.data.frame(seq(1000, 40000, 1000))
colnames(dat) <- "D"
dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
dat$YEARMARKER <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
dat$D_err <- rep(100, 40)
dat$d18Oc_err <- rep(0.1, 40)
importlist <- data_import_object(input_object = dat) # Run function on 
# attached dummy data

# Create bad data file lacking YEARMARKER column
bad_dat <- as.data.frame(seq(1000, 40000, 1000))
colnames(bad_dat) <- "D"
bad_dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
\dontrun{importlist <- data_import_object(input_object = bad_dat)}
}
