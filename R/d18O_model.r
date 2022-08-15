#' Function to convert SST data to d18O
#' 
#' Takes a matrix of SST data (in degrees C) against time (in days), information
#' about the d18O value (in permille VSMOW) of the water and how it changes
#' through the year and the transfer function used for of the record (e.g.
#' Kim and O'Neil, 1997 or Grossman and Ku, 1986). Converts the SST data to d18O
#' data using the supplied empirical transfer function.
#' 
#' @param SST Matrix with a time column (values in days) and an SST column
#' (values in degrees C)
#' @param d18Ow Either a single value (constant d18Ow) or a vector of length
#' equal to the period in SST data (365 days by default) containing information
#' about seasonality in d18Ow. Defaults to constant d18Ow of 0 permille VSMOW
#' (the modern mean ocean value)
#' @param transfer_function String containing the name of the transfer function
#' (for example: \code{"KimONeil97"} or \code{"GrossmanKu86"}). Defaults to
#' Kim and O'Neil (1997).
#' @return A vector containing d18O values for each SST value in \code{"SST"}
#' @references Grossman, E.L., Ku, T., Oxygen and carbon isotope fractionation in biogenic
#' aragonite: temperature effects, _Chemical Geology_ **1986**, _59.1_, 59-74.
#'     \doi{10.1016/0168-9622(86)90057-6}
#' Kim, S., O'Niel, J.R., Equilibrium and nonequilibrium oxygen
#' isotope effects in synthetic carbonates, _Geochimica et Cosmochimica Acta_
#' **1997**, _61.16_, 3461-3475.
#'     \doi{10.1016/S0016-7037(97)00169-5}
#' Dettman, D.L., Reische, A.K., Lohmann, K.C., Controls on the stable isotope
#' composition of seasonal growth bands in aragonitic fresh-water bivalves
#' (Unionidae), _Geochimica et Cosmochimica Acta_ **1999**, _63.7-8_, 1049-1057.
#'     \doi{10.1016/S0016-7037(99)00020-4}
#' Brand, W.A., Coplen, T.B., Vogl, J., Rosner, M., Prohaska, T., Assessment of
#' international reference materials for isotope-ratio analysis (IUPAC Technical
#' Report), _Pure and Applied Chemistry_ **2014**, _86.3_, 425-467.
#'     \doi{10.1515/pac-2013-1023}
#' Kim, S.-T., Coplen, T. B., and Horita, J.: Normalization of stable
#' isotope data for carbonate minerals: Implementation of IUPAC
#' guidelines, Geochim. Cosmochim. Ac. **2015** 158, 276-289.
#'     \doi{10.1016/j.gca.2015.02.011}
#' @examples
#' # Create dummy SST data
#' t <- seq(1, 40, 1)
#' T <- sin((2 * pi * (seq(1, 40, 1) - 8 + 10 / 4)) / 10)
#' SST <- cbind(t, T)
#' # Run d18O model function
#' d18O <- d18O_model(SST, 0, "KimONeil97")
#' @export
d18O_model <- function(SST, # Function that converts SST values into d18O 
    d18Ow = 0, # Information on the d18O of seawater (d18Ow), should be either one single number (constant) or a vector of length 365 (value for every day of the year)
    transfer_function = "KimONeil97" # Supply transfer function. Current options are: Kim and O'Neil (calcite calibration; 1997; "KimONeil97") or Grossman and Ku (aragonite calibration; 1986; "GrossmanKu86")
    ){

    if(!is.numeric(d18Ow)){
        d18Ow <- rep(0, length(SST[,1])) # Set default d18Ow vector to constant at d18Ow = 0
    }else if(length(d18Ow) == 1){
        d18Ow <- rep(d18Ow, length(SST[,1])) # Create constant d18Ow vector if only one number is given
    }else{
        d18Ow <- c(0, rep(d18Ow, length(SST[,1])/length(d18Ow))) # If d18Ow is a vector with more than one value, multiply it to reach the same length of SST (multiply by "years")
    }
    if(transfer_function == "KimONeil97"){
        d18Oc <- cbind(SST[,1], ((exp(((18.03 * 10 ^ 3) / (SST[,2] + 273.15) - 32.42) / 1000) * (((d18Ow - 30.92) / 1.03092) / 1000 + 1)) - 1) * 1000) # Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014) and Kim et al. (2015)
    }else if(transfer_function == "GrossmanKu86"){
        d18Oc <- cbind(SST[,1], (20.6 - SST[,2]) / 4.34 + d18Ow - 0.2) # Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
    }else{
        return("ERROR: Supplied transfer function is not recognized")
    }
    return(d18Oc)
}