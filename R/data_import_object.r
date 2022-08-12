#' Function to prepare d18O data from a preexisting object
#' and process yearmarkers and calculation windows
#' 
#' Takes the name of an object, reads the user-provided yearmarkers in the file
#' and uses them as a basis for the length of windows used throughout the model.
#' This ensures that windows are not too short and by default contain at least
#' one year of growth for modeling.
#' 
#' @param input_object Name of the object that contains sampling distance and
#' d18O data. Note that sampling distance should be given in micrometers,
#' because the SCEUA model underperforms when the growth rate figures are very
#' small (<0.1 mm/day).
#' @return A list containing an object with the original data and details on
#' the position and length of modeling windows
#' @examples
#' Create virtual data
#' dat <- as.data.frame(seq(1000, 40000, 1000))
#' colnames(dat) <- "D"
#' dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
#' dat$YEARMARKER <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
#' dat$D_err <- rep(100, 40)
#' dat$d18Oc_err <- rep(0.1, 40)
#' importlist <- data_import_object(input_object = dat) # Run function on 
#' # attached dummy data
#' 
#' # Create bad data file lacking YEARMARKER column
#' bad_dat <- as.data.frame(seq(1000, 40000, 1000))
#' colnames(bad_dat) <- "D"
#' bad_dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
#' \dontrun{importlist <- data_import_object(input_object = bad_dat)}
#' @export
data_import_object <- function(input_object){

    # If correct headers are included in the object, the column names don't need to be set
    # cols <- c("D", "d18Oc", "YEARMARKER", "D_err", "d18Oc_err")
    # colnames(dat) <- cols[1:length(dat[1, ])] # Give column names
    # WARNING: It is important that the columns in the datafile have the same meaning as defined here.
    # If one of the error terms (e.g. the error on the depth measurement: "D_err") is missing, it can also be added to the datafile as a column filled with zeroes (indicating the error is 0)

    # Check the structure and names of the import dataframe
    # Check if all 5 used columns are present
    if(!all(c("D", "d18Oc", "YEARMARKER", "D_err", "d18Oc_err") %in% colnames(input_object))){
        # Check if the three basic columns (without the SDs) are present
        if(!all(c("D", "d18Oc", "YEARMARKER") %in% colnames(input_object))){
            return(paste("ERROR: Input data lacks columns:", # If the three minimum requisite columns are not present, break operation
                    c("D", "d18Oc", "YEARMARKER")[
                        which(!(c("D", "d18Oc", "YEARMARKER") %in% colnames(input_object)))
                    ]
                )
            )
        }else{
            # If three basic columns are present, but error columns are missing, set errors to zero
            if(!("D_err" %in% colnames(input_object))){
                input_object$D_err <- rep(0, nrow(input_object))
                print("WARNING: D error not found, set to zero")
            }
            if(!("d18Oc_err" %in% colnames(input_object))){
                input_object$d18Oc_err <- rep(0, nrow(input_object))
                print("WARNING: d18Oc error not found, set to zero")
            }
        }
    } # No action is required if all columns are present
    
    input_object <- input_object[order(input_object$D), ] # Order data by D

    # Check for duplicate depth values and average them out
    if(TRUE %in% duplicated(input_object$D)){
        input_object <- input_object %>%
            group_by(D) %>%
            summarize(
                N = n(),
                D = median(D, na.rm = TRUE),
                d18Oc = median(d18Oc, na.rm = TRUE),
                YEARMARKER = max(YEARMARKER, na.rm = TRUE),
                d18Oc_err = sqrt(sum(d18Oc_err ^ 2, na.rm = TRUE)),
                D_err = sqrt(sum(D_err ^ 2, na.rm = TRUE)) / N
            )
        print("WARNING: Duplicated depth values were found and the median values were used")
    }

    # Define sliding window based on indicated year markers
    YEARMARKER <- which(input_object$YEARMARKER == 1) # Read position of yearmarkers in data.
    yearwindow <- diff(which(input_object$YEARMARKER == 1)) # Calculate the number of datapoints in each year between consecutive year markers
    if(length(yearwindow) > 1){
        dynwindow <- approx( # Interpolate between the numbers of annual datapoints to create list of starting positions of growth windows and their size for age modeling
            x = YEARMARKER[-length(YEARMARKER)],
            y = yearwindow,
            xout = 1:(length(input_object$D) - yearwindow[length(yearwindow)] + 1), # X indicates starting positions of windows used for age modeling
            method = "linear",
            rule = 2 # Window sizes for beginning given as NA's, for end equal to last value
        )
        dynwindow$y <- round(dynwindow$y) # Round off window sizes to integers
        dynwindow$y[dynwindow$y < 10] <- 10 # Eliminate small window sizes to lend confidence to the sinusoidal fit
        overshoot <- which((dynwindow$x + dynwindow$y) > length(input_object$D)) # Find windows that overshoot the length of data
        dynwindow$x <- dynwindow$x[-overshoot] # Remove overshooting windows
        dynwindow$y <- dynwindow$y[-overshoot] # Remove overshooting windows
        if((length(dynwindow$x) + dynwindow$y[length(dynwindow$x)] - 1) < length(input_object$D)){ # Increase length of the final window in case samples at the end are missed due to jumps in window size
            dynwindow$y[length(dynwindow$y)] <- dynwindow$y[length(dynwindow$y)] + (length(input_object$D) - (length(dynwindow$x) + dynwindow$y[length(dynwindow$x)] - 1))
        }
    }else if(length(yearwindow) == 1){ # Catch exception of datasets with only two yearmarkers
        dynwindow <- data.frame(
            x = 1:(length(input_object$D) - yearwindow[length(yearwindow)] + 1),
            y = rep(yearwindow, (length(input_object$D) - yearwindow[length(yearwindow)] + 1))
        )
    }else{
        return("ERROR: Need at least 2 year markers to estimate window size")
    }
    return(list(input_object,dynwindow))
}