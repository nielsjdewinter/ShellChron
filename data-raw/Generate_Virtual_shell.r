# Generate virtual shell data

# Function used to linearly subsample data at new depth values
subsample_mean <- function(dailydata, dailydepth, newdepth, plot=FALSE){
    newdata <- vector()
    for(i in 1:length(newdepth)){ # Loop through all new depth values (samples)
        # Find start and end positions of range to average
        if(i == 1){
            pos1 <- 1 # If first sample of the record, then start position is equal to start of the record (pos1 = 1)
            pos2 <- length(which((newdepth[i + 1] + newdepth[i]) / 2 > dailydepth)) # End position is equal to the middle between the sample depth and the depth of the next sample
        }else{
            if(i == length(newdepth)){
                pos1 <- length(which((newdepth[i - 1] + newdepth[i]) / 2 > dailydepth)) # Start position is equal to the middle between the sample depth and the depth of the previous sample
                pos2 <- length(dailydepth) # If last sample of the record, then the last position is equal to the end of the record (pos2 = length(dailydepth))
            }else{
                pos1 <- length(which((newdepth[i - 1] + newdepth[i]) / 2 > dailydepth)) # Start position is equal to the middle between the sample depth and the depth of the previous sample
                pos2 <- length(which((newdepth[i + 1] + newdepth[i]) / 2 > dailydepth)) # End position is equal to the middle between the sample depth and the depth of the next sample
            }
        }
        newdata <- append(newdata, mean(dailydata[pos1:pos2])) # Calculate the new data value for each sample by averaging the range of datapoints between start and end position
    }
    if(plot == TRUE){ # Create plot showing subsampling if requested
        dev.new()
        plot(dailydepth, dailydata, type = "l")
        points(newdepth, newdata, col = "red")
    }
    return(newdata)
}

# Function for creating d18Oc and D47 data from set of growth conditions and sampling resolution.
# All input data are vectors at daily resolution, except for D, which depends on the sampling resolution and shell length
shellmodel <- function(time, SST, GR, d18Osw, D, AV = FALSE, plot = FALSE){
    Dday <- cumsum(GR) # Create vector linking days to depth values
    if(AV == FALSE){
        SSTnew <- subsample(SST, Dday, D) # Subsample SST along the new sample set
        d18Oswnew <- subsample(d18Osw, Dday, D) # Subsample d18Osw along the new sample set
        Tynew <- subsample(Ty, Dday, D) # Subsample time (yr) along the new sample set
    }else{
        SSTnew <- subsample_mean(SST, Dday, D) # Subsample SST along the new sample set using mean values
        d18Oswnew <- subsample_mean(d18Osw, Dday, D) # Subsample d18Osw along the new sample set using mean values
        Tynew <- subsample_mean(Ty, Dday, D) # Subsample time (yr) along the new sample set using mean values
    }
    # alpha = exp((18.03*1000/(SSTnew+273.15)-33.42)/1000) # Calculate alpha of calcite fractionation
    # d18Osw_PDB = (0.97002*d18Oswnew-29.98) # Convert d18Osw to PDB
    # d18Oc = ((alpha * (d18Osw_PDB/1000 + 1)) - 1)*1000
    d18Oc <- (exp((18.03 * 1000 / (SSTnew + 273.15) - 33.42) / 1000) * ((0.97002 * d18Oswnew - 29.98) / 1000 + 1) - 1) * 1000 # Calculate d18O of calcite for each sample according to Kim and O'Neil, 1997
    D47 <- (0.0449 * 10^6) / (SSTnew + 273.15) ^ 2 + 0.167 # Calculate D47 of calcite for each sample according to Kele et al., 2015 modified by Bernasconi et al., 2018
    if(plot == TRUE){ # Create plots of new data if requested
        dev.new()
        plot(D, d18Oc, col = "blue")
        par(new = TRUE)
        plot(D, D47, axes = FALSE, bty = "n", xlab = "", ylab = "", col = "red")
        axis(side = 4, at = pretty(range(D47)))
    }
    dat <- cbind(Tynew, D, d18Oc, D47) # Combine new data for export
    return(dat) # Return the new depth, d18Oc and D47 series
}

# Set boundary conditions
Td <- seq(1, 6 * 365, 1) # Create timeline of 6 years in days
Ty <- Td / 365 # Convert to years
MAT <- 20 # Set mean annual temperature
Amp <- 10 # Set seasonal amplitude
Sext <- 2 * Amp # Calculate extent of seasonal variability
TSD <- 1.5 # Set the degree of random non-seasonal noise on the SST curve ("weather")
SST <- rnorm(length(Ty), MAT + Amp * sin(2 * pi * Ty), TSD) # Create virtual daily SST data
GR <- rep(10 / 365, length(Ty)) # Set growth rate to 10 mm/yr, create daily GR vector
DSD <- 0.6 # Set the degree of random non-seasonal noise on the d18Osw curve ("salinity fluctuations")
d18Osw <- rnorm(length(Ty), rep(0, length(Ty)), DSD) # Set d18Osw to 0 permille VSMOW, create daily d18Osw vector

SR <- 0.75 # Set uneven sampling resolutions

# Calculate virtual data
Virtual_shell <- as.data.frame(shellmodel(Ty, SST, GR, d18Osw, D, AV = TRUE))
Vritual_shell$D_err <- rep(0.1, length(Virtual_shell[, 1])) # Add uncertainties on D (in mm)
Vritual_shell$d18Oc_err <- rep(0.1, length(Virtual_shell[, 1])) # Add uncertainties on d18Oc (in permille)
# D47 data is removed and YEARMARKER column is added manually by user identification of the year transitions