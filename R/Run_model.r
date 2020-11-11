#' Function that optimizes sinusoid parameters to fit d18O data
#' 
#' The second core function of the ShellChron growth model. Loops
#' through all data windows and uses the \code{growth_model} function
#' to create d18O series that match the input data. This step is
#' iterated and optimized (minimizing the Sum of Squared Residuals)
#' through the SCEUA algorithm (by Duan et al., 1992) which finds
#' the optimal input parameters to the growth rate and Sea Surface
#' Temperature (SST) sinusoids to simulate d18O data.
#' @param dat Matrix containing the input data
#' @param dynwindow Information on the position and length of modelling
#' windows
#' @param mineral Mineralogy of record (default = "calcite")
#' @param d18Ow Either a single value (constant d18Ow) or a vector of length
#' equal to the period in SST data (365 days by default) containing information
#' about seasonality in d18Ow. Defaults to constant d18Ow of 0 permille VSMOW
#' (the modern mean ocean value)
#' @param T_per Period of SST sinusoid (in days; default = 365)
#' @param G_per Period of growth rate sinusoid (in days; default = 365)
#' @param t_int Time interval (in days; default = 1)
#' @param t_maxtemp Timing of the warmest day of the year (in julian day; 
#' default = 182.5, or May 26th halfway through the year)
#' @param MC Number of Monte Carlo simulations to apply for error propagation
#' Default = 1000
#' @param plot Should results of modelling be plotted? \code{TRUE/FALSE}
#' @return A list containing the \code{resultarray} which contains the full
#' result of all simulations on each data window and the \code{parmat} listing
#' all optimized growth rate and SST parameters used to model d18O in each data
#' window
#' @seealso Duan, Qingyun, Soroosh Sorooshian, and Vijai Gupta. "Effective and
#' efficient global optimization for conceptual rainfall runoff models." Water
#' resources research 28.4 (1992): 1015-1031. https://doi.org/10.1029/91WR02985
#' @references package dependencies: ggplot2 3.2.1; rtop 0.5.14
#' Function dependencies: sinreg, d18O_model, growth_model
#' @examples
#' # Create dummy input data column by column
#' dat <- as.data.frame(seq(1000, 40000, 1000))
#' colnames(dat) <- "D"
#' dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
#' dat$YEARMARKER <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
#' dat$D_err <- rep(100, 40)
#' dat$d18Oc_err <- rep(0.1, 40)
#' # Create dummy dynwindow data
#' dynwindow <- as.data.frame(seq(1, 29, 2))
#' colnames(dynwindow) <- "x"
#' dynwindow$y <- rep(12, 15)
#' # Run model function
#' \donttest{resultlist <- run_model(dat,
#'     dynwindow,
#'     "calcite",
#'     d18Ow = 0,
#'     T_per = 365,
#'     G_per = 365,
#'     t_int = 1,
#'     t_maxtemp = 182.5,
#'     MC = 1000,
#'     plot = FALSE)}
#' @export
run_model <- function(dat, # Core function to run the entire model on the data (dat)
    dynwindow, # The window vetor resulting from reading in the data 
    mineral = "calcite",
    d18Ow = "default",
    T_per, # Temperature sinusoid parameters
    G_per, # Growth sinusoid parameters
    t_int = 1, # Default time interval = 1 day
    t_maxtemp = 182.5, # Default time (day) at which maximum temperature is reached is 182.5 (exactly halfway through the year, or 1st of June)
    MC = 1000, # If errors = TRUE, give the number of iterations for Monte Carlo simulation used in error propagation (default = 1000, if MC = 0 no eror propagation is done)
    plot = FALSE # Should the progress of model fitting be plotted?
    ){

    d18Oc <- d18Oc_err <- Omod <- NULL # Predefine variables to circumvent global variable binding error

    # Prepare data arrays for storage of modelling results
    resultarray <- array( # Create array to contain all modelling results of overlapping windows
        rep(as.matrix(cbind(dat, matrix(NA, ncol = length(dynwindow$x), nrow = length(dat$D)))), 9), # Replicate matrix of dat + extra columns to contain all variables
        dim = c(length(dat$D), length(dynwindow$x) + length(dat[1,]), 9) # Six variables, being: Modelled d18O, residuals, Day of the Year, Growth between datapoints, Instantaneous growth rate at datapoint and Modelled temperature
    )

    parmat <- matrix(NA, nrow = 7, ncol = length(dynwindow$x)) # Matrix in which to store the modelling parameters
    colnames(parmat) <- dynwindow$x
    rownames(parmat) <- c("T_amp", "T_pha", "T_av", "G_amp", "G_pha", "G_av", "G_skw")
    
    # Prepare plot to show model progress

    if(plot == TRUE){
        dev.new()
        fitplot <- ggplot2::ggplot(dat, ggplot2::aes(D, d18Oc)) + # Create a plot showing the fit of each window on the original data, start with a plot of the original data
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = d18Oc - d18Oc_err, # Add error bars on model result (1 SD)
            ymax = d18Oc + d18Oc_err),
            width = dat$D_err) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = D - D_err,
            xmax = D + D_err),
            height = 0.05) +
        ggplot2::ggtitle("Plot of d18Oc values vs. depth. Black = data, Red = model, Errorbars = 1 SD")
        plot(fitplot)
    }

    # Estimate growth rate variability and round up to nearest multiple of 100 for conservative boundary
    GRavmax <- ceiling(max(diff(dat[dat$YEARMARKER == 1,1])) / 365 / 100) * 100

    # Find tailored range of temperatures from data
    d18Oc_range <- range(dat$d18Oc) # Find d18Oc range in data
    if(mineral == "calcite"){ # Find temperature range (to be superseded with inverse d18O_model function in later updates)
        T_range <- 18.03 * 1000 / (log((d18Oc_range - (0.97002 * rev(range(d18Ow)) - 29.98)) / 1000 + 1) * 1000 + 32.42) - 273.15 # Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014)
    }else if(mineral == "aragonite"){
        T_range <-  20.6 - 4.34 * (d18Oc_range - rev(range(d18Ow)) - 0.2) # Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
    }else{
        print("ERROR: Supplied mineralogy is not recognized")
    }
    T_max <- max(T_range)
    T_amp_max <- 2 * abs(diff(T_range))

    # Collate lower boundaries of parameters
    parl <- c(
        T_amp = 0, # Minimum T amplitude in degrees C
        T_pha = 0, # Minimum phase in days
        T_av = -4, # Minimum average T in degrees C
        G_amp = 0, # Minimum seasonal GR range in um/d
        G_pha = 0, # Minimum GR phase in days
        G_av = -1 * GRavmax, # Minimum average GR in um/d.
        G_skw = 0 # Minimum skew factor
    )

    # Collate upper boundaries of parameters
    paru <- c(
        T_amp = round(T_amp_max + 0.5, 0), # Maximum T amplitude in degrees C
        T_pha = 365, # Maximum phase in days
        T_av = round(T_max + 0.5, 0), # Maximum average T in degrees C
        G_amp = 2 * GRavmax, # Maximum seasonal GR range in um/d
        G_pha = 365, # Maximum GR phase in days
        G_av = GRavmax, # Maximum average GR in um/d (based on conservative boundaries of YEARMARKER indicators)
        G_skw = 100 # Maximum skew factor 
    )
    
    # Set parameters for SCEUA optimization
    iniflg = 1 # Flag for initial parameter array (1 = included)
    ngs = 25 # Number of complexes (sub-populations)
    maxn = 10000 # Maximum number of function evaluations allowed during optimization
    kstop = 5 # Maximum number of evolution loops before convergency
    pcento = 0.01 # Percentage change allowed in function value criterion before stop
    peps = 0.01 # Convergence level for parameter set (difference between parameters required for stop)

    # Run the model on all windows

    for(i in 1:length(dynwindow$x)){ # Loop over shell record
        print(paste("Processing Datawindow:", i, "of", length(dynwindow$x))) # Keep track of progress
        
        # Isolate year of d18O data based on window data
        Dsam <- dat[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1), 1]
        Osam <- dat[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1), 2]
        if(MC > 0){
            D_err <- dat[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1), 4] # Optional: include error on D
            O_err <- dat[dynwindow$x[i]:(dynwindow$x[i] + dynwindow$y[i] - 1), 5] # Optional: include error on d18Oc
        }else{
            D_err <- rep(0, dynwindow$y)
            O_err <- rep(0, dynwindow$y)
            MC <- 0
        }

        sinlist <- sinreg(Dsam, Osam) # Run sinusoidal regression to find initial parameter values

        # Estimate starting parameters from regression results
        O_av_start <- sinlist[[1]][1] # Export starting value for annual d18O average
        O_amp_start <- sinlist[[1]][2] # Export starting value for d18O amplitude

        if(mineral == "calcite"){
            T_av_start <- 18.03 * 1000 / (1000 * log((O_av_start - (0.97002 * mean(d18Ow) - 29.98)) / 1000 + 1) + 32.42) - 273.15  # Estimate mean temperature. Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014)
            T_amp_start <- 18.03 * 1000 / (1000 * log((O_av_start - O_amp_start - (0.97002 * mean(d18Ow) - 29.98)) / 1000 + 1) + 32.42) - 273.15 - T_av_start # Estimate temperature amplitude. Use Kim and O'Neil (1997) with conversion between VSMOW and VPDB by Brand et al. (2014)
        }else if(mineral == "aragonite"){
            T_av_start <- 20.6 - 4.34 * (O_av_start - mean(d18Ow) - 0.2) # Estimate mean temperature. Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
            T_amp_start <- 20.6 - 4.34 * (O_av_start - O_amp_start - mean(d18Ow) - 0.2) - T_av_start # Estimate mean temperature. Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
        }else{
            print("ERROR: Supplied mineralogy is not recognized")
        }

        O_pha_start <- sinlist[[1]][4] %% sinlist[[1]][3] # Estimate position (in depth of the first peak in d18O)
        O_peak <- O_pha_start + Dsam[1] # Find position of d18O peak in distance domain
        O_per_start <- sinlist[[1]][3] # Export starting value for period in distance domain
        T_pha_start <- ((O_pha_start - 0.5 * O_per_start) %% O_per_start) / O_per_start * T_per # Estimate position of first peak in temperature (low in d18O) relative to annual cycle (days)
        G_av_start <- O_per_start / G_per # Estimate average growth rate in distance/day

        years <- ceiling((diff(range(Dsam)) / O_per_start - 1) / 2) * 2 + 1 # Find next odd number to expand the number of simulated years to cover full window

        # Collate starting parameters
        par0 <- c(
            T_amp = T_amp_start,
            T_pha = T_pha_start,
            T_av = T_av_start,

            G_amp = G_av_start / 2, # Start by estimating growth rate changes by half the average
            G_pha = T_pha_start, # Start by estimating that the peak in growth rate coincides with the peak in temperature
            G_av = G_av_start,
            G_skw = 50 # Start with a no skew
        )

        invisible(capture.output( # Suppress the details on converging SCEUA
            sceua_list <- rtop::sceua(growth_model,
                par0,
                T_per = T_per,
                G_per = G_per,
                years = years,
                t_int = t_int,
                mineral = mineral,
                d18Ow = d18Ow,
                Dsam = Dsam,
                Osam = Osam,
                t_maxtemp = t_maxtemp,
                parl,
                paru,
                maxn,
                kstop,
                pcento,
                ngs,
                iniflg = iniflg,
                peps = peps,
                implicit = function(pars){sum(pars[4]/2, pars[6]) < 1} # Make sure that the cumulative GR curve is not located below 0 (if G_av < - G_amp / 2)
                )
            ))

        par1 <- sceua_list[[1]] # Eport parameters of final model
        names(par1) <- names(par0)

        result <- growth_model(par1, T_per, G_per, years, t_int, mineral, d18Ow, Dsam, Osam, t_maxtemp, plot = FALSE, MC, D_err, O_err, return = "result") # Calculate the end result of the best fit
        
        if(plot == TRUE){
            fitplot <- fitplot + # Add the new model fit to the plot to track progress of the model
                ggplot2::geom_point(data = as.data.frame(result), ggplot2::aes(Dsam, Omod), colour = "red") +
                ggplot2::geom_line(data = as.data.frame(result), ggplot2::aes(Dsam, Omod), colour = "red") 
            print(fitplot)
        }
        
        resultarray[, i + length(dat[1, ]), ] <- rbind(matrix(NA, nrow = i - 1, ncol = 9), result[, 3:11], matrix(NA, nrow = length(dat$D) - i - length(result[,1]) + 1, ncol = 9)) # Add results to result array
        parmat[, i] <- par1 # Add parameters to parameter array
    }

    # Provide names for all dimensions in the result array
    dimnames(resultarray) <- list(
        paste("sample", 1:length(resultarray[, 1, 3])),
        c(colnames(dat), paste("window", 1:length(dynwindow$x))),
        c("Modelled_d18O", "d18O_residuals", "Time_of_year", "Instantaneous_growth_rate", "Modelled temperature", "Modelled_d18O_SD", "Time_of_Year_SD", "Instantaneous_growth_rate_SD", "Modelled_temperature_SD")
    )

    colnames(parmat) <- paste("window", 1:length(parmat[1,]))
    return(list(resultarray, parmat))
}