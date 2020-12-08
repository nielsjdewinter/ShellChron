#' Function that corrects chronologies for sudden jumps in time
#' 
#' Some occurrences in the model results can lead the CumDY function
#' to detect extra year transitions, resulting in sudden jumps in
#' the shell chronology or a start of the chronology at an age
#' beyond 1 year. This function removes these sharp transitions
#' and late onset by adding or subtracting whole years to the age
#' result.
#' @param resultarray Array containing the full results of
#' the optimized growth model
#' @param T_per The period length of one year (in days)
#' @param agecorrection Correct for jumps in age (/code{TRUE}) or
#' only for starting time (/code{FALSE})
#' @param plot Should the results be plotted? (/code{TRUE/FALSE})
#' @return An updated and corrected version of \code{resultarray}
#' @references package dependencies: ggplot2 3.2.1
#' @examples
#' testarray <- array(NA, dim = c(20, 16, 9)) # Create empty array
#' # with correct third dimension
#' windowfill <- seq(10, 100, 10) # Create dummy simulation data 
#' # (ages) to copy through the array
#' for(i in 6:length(testarray[1, , 1])){
#'     testarray[, i, 3] <- c(windowfill,
#'         rep(NA, length(testarray[, 1, 3]) - length(windowfill)))
#'     windowfill <- c(NA, windowfill + 31)
#' }
#' testarray2 <- age_corr(testarray, 365, FALSE, FALSE) # Apply function on 
#' array
#' @export
age_corr <- function(resultarray,
    T_per = 365,
    plot = TRUE,
    agecorrection = TRUE
    ){
    Age <- Age_corr <- NULL # Predefine variables to circumvent global variable binding error
    # Check on glitches where consecutive windows are placed (almost) 1 year apart.
    # These are repaired by adding one T_per to the cumulative Day values of those windows in
    # resultarray[, , 3]
    mean_window_age <- data.frame(
        1:(length(resultarray[1,,1]) - 5),
        unname(colMeans(resultarray[, 6:length(resultarray[1, , 1]), 3], na.rm = TRUE)))
    colnames(mean_window_age) <- c("window", "Age")

    # Plot mean window ages
    dev.new() # Plot Window age to check for strange jumps in age
    ageplot <- ggplot2::ggplot(mean_window_age, ggplot2::aes(window, Age)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::ggtitle("Plot of average ages of modelling windows") +
        ggplot2::xlab("Window #") +
        ggplot2::scale_y_continuous("Age (days)", seq(0, 365 * ceiling(max(unname(colMeans(
            resultarray[, 6:length(resultarray[1,,1]),3], na.rm = TRUE)),
            na.rm = TRUE) / 365), 365))
    plot(ageplot)

    if(agecorrection == TRUE){
        agecorr <- rep(0, length(mean_window_age[, 2])) # Create vector to store corrections
        # (in days) to the mean ages of modelling windows
        for(i in 1:(length(mean_window_age[,2]) - 1)){ # Loop through all windows
            if(mean_window_age[i + 1, 2] - mean_window_age[i, 2] > T_per / 2){ # If there is a
            # large (> half a year) positive step between consecutive windows...
                # ...subtract one year from the subsequent windows
                agecorr[(i + 1):length(agecorr)] <- agecorr[(i + 1):length(agecorr)] - T_per
            }else if(mean_window_age[i + 1, 2] - mean_window_age[i, 2] < T_per / -2){ # If there is
                # a large (> half a year) negative drop between consecutive windows, add one year
                # from the subsequent windows
                agecorr[(i + 1):length(agecorr)] <- agecorr[(i + 1):length(agecorr)] + T_per 
            }
        }

        if(min(mean_window_age$Age + agecorr) > T_per){
            agecorr <- agecorr - 365 # Subtract whole number of years from correction if the
            # smallest value is more than one year old
        }

        mean_window_age$Age_corr <- mean_window_age$Age + agecorr
        mean_window_age$correction <- agecorr

        if(plot == TRUE){ # Plot updated result
            ageplot <- ggplot2::ggplot(mean_window_age, ggplot2::aes(window, Age)) + # Plot Window
            # age to check for strange jumps in age
            ggplot2::geom_line() +
            ggplot2::geom_point() +
            ggplot2::ggtitle("Plot of average ages of modelling windows") +
            ggplot2::xlab("Window #") +
            ggplot2::scale_y_continuous("Age (days)", seq(0, 365 * ceiling(max(unname(colMeans(
                resultarray[, 6:length(resultarray[1,,1]),3], na.rm = TRUE)),
                na.rm = TRUE) / 365), 365)) +
            ggplot2::geom_line(ggplot2::aes(window, Age_corr), col = "red") # Add corrected window
            # age to plot
            plot(ageplot)
        }
        
        resultarray[, 6:length(resultarray[1, , 1]), 3] <- resultarray[, 6:length(
            resultarray[1, , 1]), 3] + matrix(agecorr, nrow = length(resultarray[, 1, 1]),
            ncol = length(agecorr), byrow = TRUE) # Apply correction on the results of
            # age modelling
    }else{
        resultarray[, 6:length(resultarray[1, , 1]), 3] <- resultarray[, 6:length(
            resultarray[1, , 1]), 3] - floor(min(mean_window_age$Age) / 365) * 365 # Correct age if
            # all windows occur past year one
    }
    return(resultarray)
}