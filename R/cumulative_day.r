#' Function to detect year transitions and calculate cumulative age of model results
#' 
#' Takes the result of iterative growth modelling and
#' transforms data from Julian Day (0 - 365) to cumulative
#' day of the shell age by detecting where transitions
#' from one year to the next occur and adding full years
#' (365 days) to simulations in later years.
#' @param resultarray Array containing the full results of
#' the optimized growth model
#' @param plotyearmarkers Should the location of identified year
#' transitions be plotted? \code{TRUE/FALSE}
#' @param export_peakid Should the result of peak identification
#' be plotted? \code{TRUE/FALSE}
#' @param path Export path (defaults to tempdir())
#' @return A new version of the Julian Day tab of the resultarray 
#' with Julian Day model estimates replaced by estimates of 
#' cumulative age of the record in days.
#' @references package dependencies: zoo 1.8.7; scales 1.1.0; graphics
#' function dependencies: peakid
#' @importFrom graphics plot
#' @examples
#' testarray <- array(NA, dim = c(40, 36, 9)) # Create empty array
#' # with correct third dimension
#' windowfill <- seq(50, 500, 50) %% 365 # Create dummy simulation data 
#' # (ages) to copy through the array
#' for(i in 6:length(testarray[1, , 1])){
#'     testarray[, i, 3] <- c(windowfill, rep(NA, length(testarray[, 1, 3]) - 
#'         length(windowfill)))
#'     windowfill <- c(NA, (windowfill + 51) %% 365)
#' }
#' # Add dummy /code{D} column.
#' testarray[, 1, 3] <- seq(1, length(testarray[, 1, 3]), 1)
#' # Add dummy YEARMARKER column
#' testarray[, 3, 3] <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
#' # Add dummy d18Oc column
#' testarray[, 2, 3] <- sin((2 * pi * (testarray[, 1, 3] - 8 + 7 / 4)) / 7)
#' testarray2 <- suppressWarnings(cumulative_day(testarray, FALSE, FALSE, tempdir()))
#' # Apply function on array
#' @export
cumulative_day <- function(resultarray, # Align Day of year results from modelling in different windows to a common time axis
    plotyearmarkers = TRUE, # Plot peak fitting?
    export_peakid = TRUE, # Export data on how the boundaries between years were found?
    path = tempdir()
    ){
    
    dat <- resultarray[, 1:5, 3] # isolate original data
    
    # Recognition of the boundaries between years.
    # Method one: Apply sinusoidal function to the julian day simulations
    JDdat <- round(resultarray[, (length(dat[1, ]) + 1):length(resultarray[1, , 1]), 3]) # Isolate julian day simulations
    JDdat <- matrix(sin(2 * pi * (JDdat + 365/4) / 365), ncol = ncol(JDdat)) # Convert julian day to sinusoidal value (end and start of year = 1)
    JDends <- data.frame(Depth = dat[, 1],
        Yearmarker = dat[, 3],
        YEsin = scales::rescale(rowSums(JDdat, na.rm = TRUE), c(0, 1)) # Create depth series of positions which most likely represent a year end, rescale to a scale from 0 to 1
        )  
    
    # Method two: Give weight to the first and final days in reconstructions
    weightvector <- c(seq(10, 1, -1), rep(0, 346), seq(1, 10, 1)) # Create vector of weights to be given to days near the start and end of the year
    JDdat <- round(resultarray[, (length(dat[1, ]) + 1):length(resultarray[1, , 1]), 3]) # Isolate julian day simulations
    JDdat <- matrix(weightvector[JDdat + 1], ncol = ncol(JDdat)) # Apply weighting to starts and ends of the year
    JDends$YEweight <- scales::rescale(rowSums(JDdat, na.rm = TRUE), c(0, 1)) # Add normalized results to depth series

    # Method three: Use instances within the window simulations where the end of year is recorded
    JDdat <- resultarray[, (length(dat[1, ]) + 1):length(resultarray[1, , 1]), 3] # Isolate julian day simulations
    JDdat <- rbind(rep(NA, length(JDdat[1, ])), diff(JDdat) < 0) # Matrix of positions in individual windows where end of year (day 365) is recorded
    YEcount <- rowSums(JDdat, na.rm = TRUE) / max(rowSums(JDdat, na.rm = TRUE), na.rm = TRUE) # Aggregate the counted ends of years in simulations
    smoothfactor <- min(diff(which(JDends[, 2] == 1))) # Define smoothing factor for year end counts based on yearmarkers (take thickness of year with least growth to be conservative)
    YEcount <- c(rep(0, floor(smoothfactor / 2)), zoo::rollmean(YEcount, smoothfactor, align = "center"), rep(0, smoothfactor - floor(smoothfactor / 2) - 1)) # Smooth record of year end counts and pad with zeroes 
    JDends$YEcount <- scales::rescale(YEcount, c(0, 1)) # Add normalized results to depth series

    # Method four: Use maxima in d18Oc in original data
    yearpos <- c(1, which(JDends[, 2] == 1), length(JDends$Depth)) # Extract positions of yearmarkers and include start and end of record
    YE18O <- vector() # Create vector for the position of the maximum d18O value
    for(m in 1:(length(yearpos) - 1)){ # Loop through positions of yearmarkers
        if(m %in% 2:(length(yearpos) - 2)){
            maxpos <- which(dat[yearpos[m] : (yearpos[m + 1] - 1), 2] == max(dat[yearpos[m] : (yearpos[m + 1] - 1), 2])) + yearpos[m] - 1 # Find the position of the maximum value in the d18O data in that year
            if(length(maxpos) > 1){
                maxpos = round(median(maxpos)) # Prevent multiple values in maxpos (gives errors further in the calculations)
            }
            days <- seq(1, yearpos[m + 1] - yearpos[m], 1) * 365 / (yearpos[m + 1] - yearpos[m]) # Define sequence of "days" values with length = number of datapoints in the year
            sinusoid <- sin(2 * pi * (days - rep((maxpos - yearpos[m]) / (yearpos[m + 1] - yearpos[m]) * 365 - 365 / 4, length(days))) / 365) # Create sinusoid with peak at peak in d18Oc
            sinusoid[which(days > ((maxpos - yearpos[m]) / (yearpos[m + 1] - yearpos[m]) + 0.5) * 365 | days < ((maxpos - yearpos[m]) / (yearpos[m + 1] - yearpos[m]) - 0.5) * 365)] <- -1 # Assign lowest value (-1) to all datapoints more than 1/2 period away from the maximum to prevent false peaks
            YE18O <- append(YE18O, sinusoid) # add sinusoid values to running vector
        }else if(m == 1){
            maxpos <- which(dat[yearpos[m + 1] : (yearpos[m + 2] - 1), 2] == max(dat[yearpos[m + 1] : (yearpos[m + 2] - 1), 2])) + yearpos[m + 1] - 1 # Find the position of the maximum value in the d18O data of the next year (the first year that is complete)
            if(length(maxpos) > 1){
                maxpos = round(median(maxpos)) # Prevent multiple values in maxpos (gives errors further in the calculations)
            }
            days <- seq(yearpos[m] - yearpos[m + 1] + 1, yearpos[m + 2] - yearpos[m + 1], 1) * 365 / (yearpos[m + 2] - yearpos[m + 1]) # Define sequence of "days" values
            sinusoid <- sin(2 * pi * (days - rep((maxpos - yearpos[m + 1]) / (yearpos[m + 2] - yearpos[m + 1]) * 365 - 365 / 4, length(days))) / 365) # Create sinusoid with peak at peak in d18Oc
            sinusoid[which(days > ((maxpos - yearpos[m + 1]) / (yearpos[m + 2] - yearpos[m + 1]) + 0.5) * 365 | days < ((maxpos - yearpos[m + 1]) / (yearpos[m + 2] - yearpos[m + 1]) - 0.5) * 365)] <- -1 # Assign lowest value (-1) to all datapoints more than 1/2 period away from the maximum to prevent false peaks
            YE18O <- append(YE18O, sinusoid[1:(yearpos[m + 1] - yearpos[m])]) # Add sinusoid values for first bit of record to the vector
        }else if(m == (length(yearpos) - 1)){
            maxpos <- which(dat[yearpos[m - 1] : (yearpos[m] - 1), 2] == max(dat[yearpos[m - 1] : (yearpos[m] - 1), 2])) + yearpos[m - 1] - 1 # Find the position of the maximum value in the d18O data of the previous year (the last year that is complete)
            if(length(maxpos) > 1){
                maxpos = round(median(maxpos)) # Prevent multiple values in maxpos (gives errors further in the calculations)
            }
            days <- seq(1, yearpos[m + 1] - yearpos[m - 1], 1) * 365 / (yearpos[m] - yearpos[m - 1]) # Define sequence of "days" values
            sinusoid <- sin(2 * pi * (days - rep((maxpos - yearpos[m - 1]) / (yearpos[m] - yearpos[m - 1]) * 365 - 365 / 4, length(days))) / 365) # Create sinusoid with peak at peak in d18Oc
            sinusoid[which(days > ((maxpos - yearpos[m - 1]) / (yearpos[m] - yearpos[m - 1]) + 0.5) * 365 | days < ((maxpos - yearpos[m - 1]) / (yearpos[m] - yearpos[m - 1]) - 0.5) * 365)] <- -1 # Assign lowest value (-1) to all datapoints more than 1/2 period away from the maximum to prevent false peaks
            YE18O <- append(YE18O, sinusoid[(yearpos[m] - yearpos[m - 1]):(yearpos[m + 1] - yearpos[m - 1])]) # Add sinusoid values for last bit of record to the vector
        }
    }
    JDends$YE18O <- scales::rescale(YE18O, c(0, 1)) # Add normalized results to depth series
    JDends$YEcombined <- scales::rescale(rowSums(JDends[, -c(1,2)]), c(0, 1)) # Combine all four methods of peak recognition into one vector

    # Use peak ID algorhythm to find peaks in the combined vector of year end indicators to use as basis for cumulative day counting in the model results
    wmin <- round(min(diff(which(JDends[,2] == 1))) / 4) # Define starting window for peak recognition as one forth the minimum width of a year
    wmax <- round(min(diff(which(JDends[,2] == 1))) / 2) # Define maximum window for peak recognition as half the minimum width of a year
    YM <- length(which(JDends[,2] == 1)) # Extract number of years in record
    X <- list(w = vector(), p = vector()) # List for storing results
    w <- wmin # Start at minimum window size
    repeat{
        peaks <- peakid(JDends$Depth, JDends$YEcombined, w = w, span = 0.05) # Identify peaks based on current threshold
        if(length(peaks$x) == YM){ # Check if the number of years is correct and break loop if it is
            break
        }else if(w < wmax){ # Check if maximum window is reached
            X$w <- append(X$w, w) # Store window size
            X$p <- append(X$p, length(peaks$x)) # Store peak number
            w <- w + 1 # Increment window size
        }else{
            X$w <- append(X$w, w) # Store window size
            X$p <- append(X$p, length(peaks$x)) # Store peak number
            w <- max(X$w[which(abs(X$p - YM) == min(abs(X$p - YM)))]) # If no windows give the exact number of years, take the window that comes closer (prioritize larger windows in case of a tie)
            peaks <- peakid(JDends$Depth, JDends$YEcombined, w = w, span = 0.05) # Recalculate peak positions with the final chosen window size before breaking the loop
            break
        }
    }
    JDends$peakid <- rep(0, length(JDends[,1])) # Add column for peakid results
    JDends$peakid[peaks$i] <- 1 # Mark the location of peaks found in the time series

    if(plotyearmarkers == TRUE){
        dev.new()
        graphics::plot(JDends$Depth,
            JDends$YEcombined,
            type = "l",
            main = "Peak fitting results",
            xlab = "Depth",
            ylab = "Probability of end of year") # Plot aggregate of yearmarkers
        points(peaks$x, rep(1, length(peaks$x)), col = "red") # Plot location of yearmarkers
    }

    if(export_peakid == TRUE){
        write.csv(JDends, file.path(path, "peakid.csv"))
    }

    # Apply calculation of years in the model on the simulation results
    JDdat <- resultarray[, (length(dat[1, ]) + 1):length(resultarray[1, , 1]), 3] # Isolate julian day simulations
    for(col in 1:length(JDdat[1, ])){ # Loop through all simulation windows
        window <- JDdat[, col] # Isolate simulation window
        peakx <- head(which(peaks$x %in% dat[which(!is.na(window)), 1]), 1) # Find position of the first year end in the window column
        if(length(peakx) == 0){
            if(peaks$x[1] < dat[head(which(!is.na(window)), 1), 1]){ # If no global year transitions fall within the window, check if there are global year transitions before the window and find the last year transition before the first value in the window
                peakx <- max(which(peaks$x < dat[head(which(!is.na(window)), 1), 1])) 
                window <- window + peakx * 365 # And add the number of years associated with this last transition to all window values
            } # If all year transitions happen after the window, no year number needs to be added to the window values
        }else{
            JDpeak <- JDdat[tail(which(dat[, 1] == peaks$x[peakx] & !is.na(window)), 1), col] # Find JD simulation belonging to the first year transition
            if(JDpeak < 182.5){ 
                window <- window + peakx * 365 # If year transition co-occurs with first half of the simulated year, all simulated values are assumed to belong to the next year
            }else{
                window <- window + (peakx - 1) * 365 # If year transition co-occurs with last half of the simulated year, all simulated values are assumed to belong to the previous year
            }
        }
        if(length(which(diff(window) < 0)) > 0){ # Check if year transitions occur within the simulation
            for(i in which(diff(window) < 0)){
                if(length(peakx) > 0){ # Check if global year transitions occur within the window
                    if(i < which(dat[, 1] == peaks$x[peakx])){
                        window[1:i] <- window[1:i] - 365 # If the transition happens before the reference point (peakx), subtract a year's worth of days from the values before to prevent adding a year twice.
                    }
                }else{
                    window[(i + 1):length(window)] <- window[(i + 1):length(window)] + 365 # Add one year's worth of days to simulations after each transition
                }
            }
        }
        JDdat[, col] <- window # Update the julian day data with the new cumulative day simulations
    }

    result <- resultarray[, , 3]
    result[, 6:length(result[1, ])] <- JDdat # Add the updated cumulative julian day results to the resultarray and export
    return(result)
}