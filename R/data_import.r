#' Function to import d18O data and process yearmarkers and calculation windows
#' 
#' Takes the name of a file that is formatted according to the standard format
#' and converts it to an object to be used later in the model. In doing so, the
#' function also reads the user-provided yearmarkers in the file and uses them
#' as a basis for the length of windows used throughout the model. This ensures
#' that windows are not too short and by default contain at least one year of
#' growth for modelling.
#' 
#' @param file_name Name of the file that contains d18O data
#' @return A list containing an object with the original data and details on
#' the position and length of modelling windows
#' @examples
#' importlist <- data_import(file_name = system.file("extdata",
#'     "Virtual_shell.csv", package = "ShellChron")) # Run function on attached
#'     # dummy data
#' @export
data_import <- function(file_name){
    dat <- read.csv(file_name, header = T) # Read in data file
    
    # If correct headers are included in the file, the column names don't need to be set
    # cols <- c("D", "d18Oc", "YEARMARKER", "D_err", "d18Oc_err")
    # colnames(dat) <- cols[1:length(dat[1, ])] # Give column names
    # WARNING: It is important that the columns in the datafile have the same meaning as defined here.
    # If one of the error terms (e.g. the error on the depth measurement: "D_err") is missing, it can also be added to the datafile as a column filled with zeroes (indicating the error is 0)
    
    dat <- dat[order(dat[, 1]),] # Order data by D

    # Define sliding window based on indicated year markers
    YEARMARKER <- which(dat$YEARMARKER == 1) # Read position of yearmarkers in data.
    yearwindow <- diff(which(dat$YEARMARKER == 1)) # Calculate the number of datapoints in each year between consecutive year markers
    dynwindow <- approx( # Interpolate between the numbers of annual datapoints to create list of starting positions of growth windows and their size for age modelling
        x = YEARMARKER[-length(YEARMARKER)],
        y = yearwindow,
        xout = 1:(length(dat$D) - yearwindow[length(yearwindow)] + 1), # X indicates starting positions of windows used for age modelling
        method = "linear",
        rule = 2 # Window sizes for beginning given as NA's, for end equal to last value
    )
    dynwindow$y <- round(dynwindow$y) # Round off window sizes to integers
    dynwindow$y[dynwindow$y < 10] <- 10 # Eliminate small window sizes to lend confidence to the sinusoidal fit
    overshoot<-which(dynwindow$x + dynwindow$y > length(dat[,1])) # Find windows that overshoot the length of dat
    dynwindow$x <- dynwindow$x[-overshoot] # Remove overshooting windows
    dynwindow$y <- dynwindow$y[-overshoot] # Remove overshooting windows
    if((length(dynwindow$x) + dynwindow$y[length(dynwindow$x)] - 1) < length(dat[, 1])){ # Increase length of the final window in case samples at the end are missed due to jumps in window size
        dynwindow$y[length(dynwindow$y)] <- dynwindow$y[length(dynwindow$y)] + (length(dat[, 1]) - (length(dynwindow$x) + dynwindow$y[length(dynwindow$x)] - 1))
    }
    return(list(dat,dynwindow))
}