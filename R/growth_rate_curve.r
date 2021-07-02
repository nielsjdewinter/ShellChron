#' Function that creates a skewed sinusoidal growth rate (GR) curve
#' from a list of parameters
#' 
#' Takes the specified parameters for amplitude, period, phase, average value
#' and skewness factor as well as the number of years specified and the time 
#' interval. It then creates a skewed sinusoid based on the boundary conditions.
#' The skewness factor (\code{G_skw}) determines whether the sinusoid is skewed
#' towards the front (\code{G_skw < 50}) or the back of the annual peak in
#' growth rate (\code{G_skw > 50}).
#' Used as intermediate step during iterative modeling.
#' 
#' @param G_par List of four parameters describing (in order) amplitude
#' (\code{G_amp}; in micrometer/day), period (\code{G_per}; in days), phase
#' (\code{G_pha} in day of the year), average growth rate (\code{G_av}; in 
#' micrometer/day) and the skewness factor (\code{G_skw} between 0 and 100)
#' @param years Length of the preferred sinusoid in number of years (defaults
#' to 1)
#' @param t_int Time interval of sinusoidal record (in days)
#' @return A matrix containing columns for time (in days) and GR
#' (in micrometer/day)
#' @references
#'   \doi{10.1016/j.palaeo.2017.09.034}
#' @examples
#' # Set parameters
#' G_amp <- 20
#' G_per <- 365
#' G_pha <- 100
#' G_av <- 15
#' G_skw <- 70
#' G_par <- c(G_amp, G_per, G_pha, G_av, G_skw)
#' # Run GR model function
#' GR <- growth_rate_curve(G_par, 1, 1)
#' @export
growth_rate_curve <- function(G_par, # Function to create growth rate as function of time based on parameters
    years = 1, # Default number of years of the record = 1
    t_int = 1 # Default time interval = 1 day
    ){

    G_amp <- G_par[1] # Seasonal range in Growth rate (um/d)
    G_per <- G_par[2] # Period of growth rate seasonality (days)
    G_pha <- G_par[3] # Timing of peak in growth rate (day of the year)
    G_av <- G_par[4] # Annual average growth rate (um/d)
    G_skw <- G_par[5] # Skewness factor in growth rate sinusoid (-)

    t <- seq(0, G_per, t_int) # Define time axis based on the number of days in a year and the time interval

    GR <- rep(NA, length(t)) # Create empty growth rate vector
    # Build GR vector piece by piece based on parameters
    # Check if t is between the time of maximum growth and the subsequent time of minimum growth but before the peak, add one period's length if it is
    if(length(t[(((t - G_pha) %% G_per) < (G_per * (100 - G_skw) / 100)) & (t < G_pha)]) > 0){ # Catch "replacement has length zero" errors
        GR[(((t - G_pha) %% G_per) < (G_per * (100 - G_skw) / 100)) & (t < G_pha)] <- G_av + G_amp/2 * sin(2 * pi * (t[(((t - G_pha) %% G_per) < (G_per * (100 - G_skw) / 100)) & (t < G_pha)] + G_per - G_pha + (G_per * (100 - G_skw) / 50)/4) / (G_per * (100 - G_skw) / 50))
    }
    # Check if t is between the time of maximum growth and the subsequent time of minimum growth but still after the peak
    if(length(t[(((t - G_pha) %% G_per) < (G_per * (100 - G_skw) / 100)) & (t >= G_pha)]) > 0){ # Catch "replacement has length zero" errors
        GR[(((t - G_pha) %% G_per) < (G_per * (100 - G_skw) / 100)) & (t >= G_pha)] <- G_av + G_amp/2 * sin(2 * pi * (t[(((t - G_pha) %% G_per) < (G_per * (100 - G_skw) / 100)) & (t >= G_pha)] - G_pha + (G_per * (100 - G_skw) / 50)/4) / (G_per * (100 - G_skw) / 50))
    }
    # Check if t is between the time of maximum growth and the previous time of minimum growth but still before the peak, subtract one period's length if it is
    if(length(t[(((t - G_pha) %% G_per) >= (G_per * (100 - G_skw) / 100)) & (t > G_pha)]) > 0){
        GR[(((t - G_pha) %% G_per) >= (G_per * (100 - G_skw) / 100)) & (t > G_pha)] <- G_av + G_amp/2 * sin(2 * pi * (t[(((t - G_pha) %% G_per) >= (G_per * (100 - G_skw) / 100)) & (t > G_pha)] - G_per - G_pha + (G_per * G_skw / 50)/4) / (G_per * G_skw / 50))
    }
    if(length(t[(((t - G_pha) %% G_per) >= (G_per * (100 - G_skw) / 100)) & (t <= G_pha)]) > 0){
        GR[(((t - G_pha) %% G_per) >= (G_per * (100 - G_skw) / 100)) & (t <= G_pha)] <- G_av + G_amp/2 * sin(2 * pi * (t[(((t - G_pha) %% G_per) >= (G_per * (100 - G_skw) / 100)) & (t <= G_pha)] - G_pha + (G_per * G_skw / 50)/4) / (G_per * G_skw / 50))
    }
    # Replace negative growth rates with zeroes
    GR[GR < 0] <- 0
    # Multiply t and GR with number of years
    if(years > 1){
        t <- c(t, rep(t[-1], years - 1) + rep(seq(365, 365 * (years - 1), 365), each = 365))
        GR <- c(GR, rep(GR[-1], years - 1))
    }
    # Collate results and export
    res <- cbind(t, GR)
    return(res)
}