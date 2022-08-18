
#' Function BLATANTLY STOLEN FROM SHELLCHRON PACKAGE to convert d18O data to SST
#'
#' Takes a matrix of d18Os data (in permille VPDB) against distance measures (in any unit), information
#' about the d18O value (in permille VSMOW) of the water and how it changes
#' from one value to another, and the transfer function used for of the record (e.g.
#' Kim and O'Neil, 1997 or Grossman and Ku, 1986). Converts the d18O data to SST
#' data (in degrees Celsius) using the supplied empirical transfer function.
#'
#' @param d18Oc Matrix with a distance column (values in any unit) and an d18Oc column
#' (values in permille VPDB)
#' @param T Temperature of seawater in degrees C. Either a single value (constant temperature)
#' or a vector of length equal to the number of d18O values
#' @param transfer_function String containing the name of the transfer function
#' (for example: \code{"KimONeil97"} or \code{"GrossmanKu86"}). Defaults to
#' Kim and O'Neil (1997).
#' @return A vector containing SST values for each d18Os value in \code{"d18Oc"}
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
#' @examples
#' # Create dummy d18Os data
#' dist <- seq(1, 40, 1) #distance
#' val <- sin((0.5 * pi * (dist)) / 5)+1 # d18Os
#' d18O <- cbind(dist, val)
#' # Run SST model function
#' SST <- SST_model(d18O, 0, "KimONeil97")
#' @export
d18Ow_from_d18O_T <- function(d18Oc, # Function that calculates d18Ow from d18Oc and temperature values
                       T, # Information on the temperature of seawater in degrees C, should be either one single number (constant) or a vector of length identical to the length of d18Os values.
                       transfer_function = "KimONeil97" # Supply transfer function. Current options are: Kim and O'Neil (calcite calibration; 1997; "KimONeil97") or Grossman and Ku (aragonite calibration; 1986; "GrossmanKu86")
){

  if(!is.numeric(T)){
    T <- rep(0, length(d18Oc[, 1])) # Set default T vector to constant at T = 0
  }else if(length(T) == 1){
    T <- rep(T, length(d18Oc[, 1])) # Create constant T vector if only one number is given
  }else{
    T <- c(0, rep(T, length(d18Oc[, 1]) / length(T))) # If T is a vector with more than one value, multiply it to reach the same length of SST (multiply by "years")
  }
  if(transfer_function == "KimONeil97"){
    d18Ow <- cbind(d18Oc[, 1], ((d18Oc[, 2] / 1000 + 1) / exp(((18.03 * 10 ^ 3) / (T - 273.15) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92)
  }else if(transfer_function == "GrossmanKu86"){
    d18Ow <- cbind(d18Oc[, 1], d18Oc[, 2] - (20.6 - T) / 4.34 + 0.2) # Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
  }else{
    return("ERROR: Supplied transfer function is not recognized")
  }
  return(SST)
}
