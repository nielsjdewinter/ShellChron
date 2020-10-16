#' Function to convert SST data to d18O
#' 
#' Takes a matrix of SST data (in degrees C) against time (in days), information
#' about the d18O value (in permille VSMOW) of the water and how it changes
#' through the year and the mineralogy of the record (calcite or aragonite). 
#' Converts the SST data to d18O data using an empirical transfer function
#' (in function of the mineralogy)
#' 
#' @param SST Matrix with a time column (values in days) and an SST column
#' (values in degrees C)
#' @param d18Ow Either a single value (constant d18Ow) or a vector of length
#' equal to the period in SST data (365 days by default) containing information
#' about seasonality in d18Ow. Defaults to constant d18Ow of 0 permille VSMOW
#' (the modern mean ocean value)
#' @param mineral String containing the name of the mineralogy (either 
#' \code{"calcite"} or \code{"aragonite"}). Defaults to calcite.
#' @return A vector containing d18O values for each SST value in \code{"SST"}
#' @examples
#' # Create dummy SST data
#' t <- seq(1, 40, 1)
#' T <- sin((2 * pi * (seq(1, 40, 1) - 8 + 10 / 4)) / 10)
#' SST <- cbind(t, T)
#' # Run d18O model function
#' d18O <- d18O_model(SST, 0, "calcite")
#' @export
d18O_model <- function(SST, # Function that converts SST values into d18O 
    d18Ow = 0, # Information on the d18O of seawater (d18Ow), should be either
    # one single number (constant) or a vector of length 365 (value for every
    # day of the year)
    mineral = "calcite" # Using either Kim and O'Neil (calcite calibration;
    # 1997) or Grossman and Ku (aragonite calibration; 1986) 
    ){

    if(!is.numeric(d18Ow)){
        d18Ow <- rep(0, length(SST[,1])) # Set default d18Ow vector to constant
        # at d18Ow = 0
    }else if(length(d18Ow) == 1){
        d18Ow <- rep(d18Ow, length(SST[,1])) # Create constant d18Ow vector if
        # only one number is given
    }else{
        d18Ow <- c(0, rep(d18Ow, length(SST[,1])/length(d18Ow))) # If d18Ow is
        # a vector with more than one value, multiply it to reach the same
        # length of SST (multiply by "years")
    }
    if(mineral == "calcite"){
        d18Oc <- cbind(SST[,1], (exp((18.03 * 1000 / (SST[,2] + 273.15) -
            32.42) / 1000) - 1) * 1000 + (0.97002 * d18Ow - 29.98)) # Use
            # Kim and O'Neil (1997) with conversion between VSMOW and VPDB by
            # Brand et al. (2014)
    }else if(mineral == "aragonite"){
        d18Oc <- cbind(SST[,1], (20.6 - SST[,2]) / 4.34 + d18Ow + 0.2) # Use
        # Grossmann and Ku (1986) modified by Dettmann et al. (1999)
    }else{
        print("ERROR: Supplied mineralogy is not recognized")
    }
    return(d18Oc)
}