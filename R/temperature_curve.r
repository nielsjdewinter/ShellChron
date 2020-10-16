#' Function that creates a sinusoidal Sea Surface Temperature (SST) curve
#' from a list of parameters
#' 
#' Takes the specified parameters for amplitude, period, phase and average value
#' as well as the number of years specified and the time interval. It then 
#' creates a sinusoid based on the boundary conditions. Used as intermediate
#' step during iterative modelling.
#' 
#' @param T_par List of four parameters describing (in order) amplitude
#' (\code{T_amp}; in degrees C), period (\code{T_per}; in days), phase
#' (\code{T_pha} in day of the year) and average temperature (\code{T_av};
#' in degrees C)
#' @param years Length of the preferred sinusoid in number of years (defaults
#' to 1)
#' @param t_int Time interval of sinusoidal record (in days)
#' @return A matrix containing columns for time (in days) and SST (in degrees C)
#' @examples
#' # Set parameters
#' T_amp <- 20
#' T_per <- 365
#' T_pha <- 150
#' T_av <- 15
#' T_par <- c(T_amp, T_per, T_pha, T_av)
#' SST <- temperature_curve(T_par, 1, 1) # Run the function
#' @export
temperature_curve <- function(T_par, # Function to create temperature as function of time based on sinusoidal parameters
    years = 1, # Default number of years of the record = 1
    t_int = 1 #  Default time interval = 1 day
    ){

    T_amp <- T_par[1] # Seasonal SST range (degr. C)
    T_per <- T_par[2] # Period of SST seasonality (days)
    T_pha <- T_par[3] # Timing of peak in SST sinusoid (day of the year)
    T_av <- T_par[4] # Annual average SST (degr. C)

    t <- seq(0, years * T_per, t_int) # Define time axis based on the number of days in a year and the time interval. Model multiple years (default = 3) to accommodate changes in growth rate over lifetime.

    SST <- T_av + T_amp/2 * sin((2 * pi * (t - T_pha + T_per/4)) / T_per) # Create SST sinusoid based on parameters and time vector

    # Collate results and export
    res <- cbind(t, SST)
    return(res)
}