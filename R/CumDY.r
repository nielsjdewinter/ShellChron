#' Function to detect year transitions and calculate cumulative age
#' of model results
#' 
#' Takes the result of iterative growth modelling and
#' transforms data from Julian Day (0 - 365) to cumulative
#' day of the shell age by detecting where transitions
#' from one year to the next occur and adding full years
#' (365 days) to simulations in later years.
#' @param resultarray Array containing the full results of
#' the optimized growth model
#' @param threshold Artificial threshold value used to recognize
#' peaks in occurrences of year transitions (default = 5)
#' @param plotyearmarkers Should the location of identified year
#' transitions be plotted? \code{TRUE/FALSE}
#' @return A new version of the resultarray with Julian Day model 
#' estimates replaced by estimates of cumulative age of the record
#' in days.
#' @references package dependencies: zoo 1.8.7
#' @examples
#' testarray <- array(NA, dim = c(20, 16, 9)) # Create empty array
#' # with correct third dimension
#' windowfill <- seq(50, 500, 50) # Create dummy simulation data 
#' # (ages) to copy through the array
#' for(i in 6:length(testarray[1, , 1])){
#'     testarray[, i, 3] <- c(windowfill, rep(NA, length(testarray[, 1, 3]) -
#'         length(windowfill)))
#'     windowfill <- c(NA, (windowfill + 51) %% 365)
#' }
#' testarray[, 1, 3] <- seq(1, length(testarray[, 1, 3]), 1) # Add
#' # dummy /code{D} column.
#' testarray2 <- cumdy(testarray, 3, FALSE) # Apply function on array
#' @export
cumdy <- function(resultarray, # Align Day of year results from modelling in different windows to a common time axis
    threshold = 5, # Threshold for separating peaks in year changes for marking the transitions between years
    plotyearmarkers
    ){ 
    
    dat <- resultarray[, 1:5, 3]
    Yearends <- rbind(rep(NA, length(resultarray[1, -(1:length(dat[1, ])), 3])), diff(resultarray[, -(1:length(dat[1, ])), 3]) < 0) # Matrix of positions in individual windows where end of year (day 365) is recorded
    Tyearmarkers <- cbind(c(0, dat[-nrow(dat), 1]) + c(0, diff(dat[, 1])), rowSums(Yearends, na.rm = TRUE)) # Aggregate of end of year (Day = 365) positions across windows against depth of record
    Tyearmarkers <- cbind(Tyearmarkers, c(rep(0, floor(threshold / 2)), zoo::rollmean(Tyearmarkers[, 2], threshold, align = "center"), rep(0, threshold - floor(threshold / 2) - 1))) # Add moving average of threshold value, pad with zeroes to match column length
    if(!all(Tyearmarkers[, 3] < (3 / threshold))){Tyearmarkers[which(Tyearmarkers[, 3] < (3 / threshold)), 3] <- 0} # Remove small numbers consisting of averages of 1 or 2 yearends in windows except in cases with very low resolution

    pks <- which(diff(sign(diff(Tyearmarkers[, 3], na.pad = FALSE)), na.pad = FALSE) < 0) + 2 # Find peaks in the number of yearmarkers by taking the second derivative
    Tyearmarkers <- cbind(Tyearmarkers, rep(0, length(Tyearmarkers[, 3]))) # Add column to store yearmarkers based on age modelling
    for(i in 1:(length(which(diff(pks) > threshold)) + 1)){ # Loop through the instances where peaks are far enough apart to be taken as separate (as judged through the threshold)
        Tyearmarkers[mean(pks[(c(0, which(diff(pks) > threshold), length(pks))[i] + 1) : c(0, which(diff(pks) > threshold), length(pks))[i + 1]]), 4] <- 1 # Combine peaks that cluster together and add markers to the mean positions where windows record the end of the year
    }
    Tyearends <- cbind(seq(1, sum(Tyearmarkers[, 4]), 1), Tyearmarkers[which(Tyearmarkers[, 4] == 1), 1]) # Find numbers and positions of mean yearmarkers in record based on age modelling

    if(plotyearmarkers == TRUE){
        dev.new()
        plot(Tyearmarkers[, c(1, 3)], type = "l") # Plot aggregate of yearmarkers
        points(Tyearmarkers[which(Tyearmarkers[, 4] == 1), c(1, 3)], col = "red") # Plot location of yearmarkers
    }

    Yearends[which(Yearends == TRUE)] <- apply(abs(outer(Tyearends[, 2], Tyearmarkers[which(Yearends == TRUE) %% nrow(Yearends), 1], "-")), 2, which.min)  # Find distance values for year ends in all windows and replace the positions of the year ends with the number of years along the records based on the nearest peak in yearends found earlier

    for(col in 1:ncol(Yearends)){ # Loop through columns and create matrix of the cumulative year in which the datapoints of all models are set.
        Yearends[min(which(!is.na(Yearends[, col]))) - 1, col] <- 0 # Replace the last NA before the modelled values start with a zero to compensate for the missing first value due to diff() function above (first line of function)
        if(length(Yearends[which(Yearends[, col] > 0), col]) > 0){ # Check if there are year ends in the column
            X <- Yearends[which(Yearends[, col] > 0), col] # Find values associated with the year ends
            row <- which(Yearends[, col] %in% X) # Find rows in which these year end values are contained
            row0 <- which(Yearends[, col] == 0) # Find rownumbers of zeroes in column
            if(length(row) > 1){ # Check if multiple year ends are present on the column
                rowcomp <- outer(row0, row, "-") # Match rownumbers of all zeroes in the column with the rownumbers of all year end values
                rowcomp[rowcomp < 0] <- NA # Remove all instances where zeroes precede year ends
                Yearends[row0, col] <- as.numeric(apply(rowcomp, 1, which.min)) # Replace zeroes with closest year end directly above
                Yearends[row0[which(is.na(Yearends[row0, col]))], col] <- X[1] - 1 # Replace all zeroes above the year end with the number of the previous year (year end - 1)
            }else{ # Alternatively, if there is only one year end:
                Yearends[row0[which(row0 > row)], col] <- Yearends[row, col] # Replace all zeroes below the year end with the number belonging to the year end
                Yearends[row0[which(row0 < row)], col] <- Yearends[row, col] - 1 # Replace all zeroes above the year end with the number of the previous year (year end - 1)
            }
        }else if(col > 1){
            Yearends[which(!is.na(Yearends[, col])), col] <- rep(as.numeric(names(sort(table(Yearends[, (col - 1)]), decreasing = TRUE)[1])), length(which(!is.na(Yearends[, col])))) # If no yearends are present, replace all zeroes with the most common year in the previous column
        }
    }

    result <- resultarray[, , 3] + cbind(matrix(0, ncol = length(dat[1, ]), nrow = length(Yearends[,1])), Yearends) * 365 # Create matrix of cumulative days in model centered on position of highest summer value (default t_maxtemp = day 182.5) and based on length of day (default T_per = G_per = 365)
    return(result)
}