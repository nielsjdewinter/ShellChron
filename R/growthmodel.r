#' Function that models a d18O curve through SST and GR sinusoids
#' 
#' The core function of the ShellChron growth model. Uses growth
#' rate and SST (Sea Surface Temperature) sinusoids to model d18O
#' data to be matched with the input. In the ShellChron modelling 
#' routine, this function is optimized using the SCEUA algorithm
#' and applied on sliding windows through the dataset to estimate
#' the age of each datapoint
#' @param pars List of parameters for temperature and growth rate sinusoids
#' \code{pars <- c(T_amp, T_pha, T_av, G_amp, G_pha, G_av, G_skw)}
#' @param T_per Period of SST sinusoid (in days; default = 365)
#' @param G_per Period of growth rate sinusoid (in days; default = 365)
#' @param years Number of years to be modelled (default = 1)
#' @param t_int Time interval (in days; default = 1)
#' @param mineral Mineralogy of record (default = "calcite")
#' @param d18Ow Either a single value (constant d18Ow) or a vector of length
#' equal to the period in SST data (365 days by default) containing information
#' about seasonality in d18Ow. Defaults to constant d18Ow of 0 permille VSMOW
#' (the modern mean ocean value)
#' @param Dsam Vector of \code{D} values serving as input (keep unit consistent
#' throughout model)
#' @param Osam Vector of \code{d18Oc} values serving as input (in permille VPDB)
#' @param t_maxtemp Timing of the warmest day of the year (in julian day; 
#' default = 182.5, or May 26th halfway through the year)
#' @param plot Should results of modelling be plotted? \code{TRUE/FALSE}
#' @param MC Number of Monte Carlo simulations to apply for error propagation
#' Default = 1000
#' @param D_err OPTIONAL: Vector containing errors on \code{Dsam}
#' @param O_err OPTIONAL: Vector containing errors on \code{Osam}
#' @param return String indicating whether to return just the Sum of Squared
#' Residuals ("SSR") or a matrix containing the results of the model and the
#' propagated uncertainties (if applicable)
#' @return Depending on the value of the "return" parameter either a single
#' value representing the Sum of Squared Residuals ("SSR") as a measure for
#' the closeness of the match between modelled d18O and input values, or a
#' matrix containing the full result of the modelling including propagated
#' uncertainties if applicable.
#' @references package dependencies: ggplot2 3.2.1
#' function dependencies: SSTcurve, d18Omodel, GRcurve, MCerr_orth
#' @examples
#' # Set parameters
#' G_amp <- 20
#' G_per <- 365
#' G_pha <- 100
#' G_av <- 15
#' G_skw <- 70
#' T_amp <- 20
#' T_per <- 365
#' T_pha <- 150
#' T_av <- 15
#' pars <- c(T_amp, T_pha, T_av, G_amp, G_pha, G_av, G_skw)
#' d18Ow <- 0
#' # Create dummy data
#' Dsam <- seq(1, 40, 1)
#' Osam <- sin((2 * pi * (seq(1, 40, 1) - 8 + 30 / 4)) / 30)
#' # Test returning residual sum of squares for optimization
#' SSR <- growthmodel(pars, T_per, G_per, Dsam = Dsam, Osam = Osam, return = "SSR")
#' # Test returning full model result
#' resmat <- growthmodel(pars, T_per, G_per, Dsam = Dsam, Osam = Osam, return = "result")
#' @export
growthmodel <- function(pars, # Growth model function to optimize using sceua
    T_per = 365, # Period of SST sinusoid (default = 365 days)
    G_per = 365, # Period of growth sinusoid (default = 365 days)
    years = 1, # Default year = 1
    t_int = 1, # Time interval, default = 1 day
    mineral = "calcite", # Mineralogy of the material (calcite and aragonite supported, default = calcite)
    d18Ow = "default", # d18Ow vector, default = constant at 0 permille VSMOW
    Dsam, # Depth data
    Osam, # d18O data
    t_maxtemp = 182.5, # Time (day) at which maximum temperature is reached. Default is 182.5 (exactly halfway through the year, or 1st of June)
    plot = FALSE, # Plot results?
    MC = 1000, # Number of MC simulations to include measurement error into error analysis. Default = 1000 (error on D and d18O measurements not considered)
    D_err = NULL, # Optional: include standard deviation of uncertainty on depth measurements, default = empty (no error incorporated in the calculation, same effect as MC = 0)
    O_err = NULL, # Optional: include standard deviation of uncertainty on d18O measurements, default = empty (no error incorporated in the calculation, same effect as MC = 0)
    return = "SSR" # What to return? Default = the total sum of squares to be minimized by the SCEUA algorithm
    ){

    # Extract variables
    T_amp <- pars[1]
    T_pha = pars[2]
    T_av = pars[3]

    G_amp = pars[4]
    G_pha = pars[5]
    G_av = pars[6]
    G_skw = pars[7]

    T_par <- c(T_amp, T_per, T_pha, T_av)
    G_par <- c(G_amp, G_per, G_pha, G_av, G_skw)

    # Run functions for SST and GR
    SST <- SSTcurve(T_par, years, t_int)
    d18Oc <- d18Omodel(SST, d18Ow, mineral)
    GR <- GRcurve(G_par, years, t_int)

    # Create distance (D) vector from GR
    D <- cumsum(GR[,2]) + Dsam[1]

    # Compare model result with data
    Omod <- approx( # Interpolate modelled d18Oc values to positions along the record.
        x = D,
        y = d18Oc[,2],
        xout = Dsam,
        method = "linear",
        rule = 2
    )

    residuals <- Osam - Omod$y # Calculate residuals

    if(return == "SSR"){
        SSR <- sum(residuals^2) # Calculate sum of squared residuals
        return(SSR)
    }else{
        # Compare model result with data
        t <- approx( # Interpolate modelled time values to positions along the record.
            x = D,
            y = d18Oc[,1],
            xout = Dsam,
            method = "linear",
            rule = 2
        )

        gr <- approx( # Interpolate modelled growth rate values to positions along the record.
            x = D,
            y = GR[,2],
            xout = Dsam,
            method = "linear",
            rule = 2
        )

        Tmod <- approx( # Interpolate modelled temperature values to positions along the record.
            x = D,
            y = SST[,2],
            xout = Dsam,
            method = "linear",
            rule = 2
        )

        TY <- (t$y - T_pha + t_maxtemp) %% 365 # Calculate time of year relative to the first T maximum in the record (T_pha) and the assumed day of the year when it occurs (default = 182.5)

        # Prepare export matrix
        resmat <- cbind(Dsam, Osam, Omod$y, residuals, TY, gr$y, Tmod$y) # Calculate matrix of end results
        colnames(resmat) <- c("Dsam", "Osam", "Omod", "Residuals", "Time_of_year", "Instantaneous_growth_rate", "Modelled_temperature")
        
        # Optional: Incorporate measurement errors on D and d18O
        if(MC > 0){
            print("Propagating measurement uncertainties")
            if(is.null(D_err)){ # If no error on D is given, but error analysis is asked (MC > 0), then D_err is zero for all samples
                D_err <- rep(0, length(Dsam))
            }
            if(is.null(O_err)){ # If no error on d18Oc is given, but error analysis is asked (MC > 0), then O_err is zero for all samples
                O_err = rep(0, length(Osam))
            }
            # Propagate the combined effect of error on D and d18Oc on the model fit
            D_err_comb <- MCerr_orth(Dsam, D_err, Osam, O_err, D, d18Oc, MC) # Combine error on D and d18Oc on the depth domain through orthogonal projection on the modelled D-d18Oc curve

            Drange <- cbind((Dsam - D_err_comb) %% D[length(D)], (Dsam + D_err_comb) %% D[length(D)]) # Find upper and lower boundaries of D error (1 SD)
            Drange_pos <- apply(abs(outer(Drange, D, FUN = "-")), c(1, 2), which.min) # Find positions of D ranges

            # Approximate SDs of d18Oc, t, gr and Tmod from range of D values within +/- 1 SD
            Omod_SD <- vector(length = length(Drange_pos[,1]))
            t_SD <- vector(length = length(Drange_pos[,1]))
            gr_SD <- vector(length = length(Drange_pos[,1]))
            Tmod_SD <- vector(length = length(Drange_pos[,1]))
            for(j in 1:length(Drange_pos[,1])){ # Loop through all positions in the window and propagate the D error onto the results
                if(Drange_pos[j, 2] < Drange_pos[j, 1]){ # If the D range contains the year edge, incorporate the outer edges of the year.
                # Find the range of modelled values contained in the 1 SD range of D error and approximate SD of modelled values from in and max values
                    Omod_SD[j] <- (max(d18Oc[c(Drange_pos[j, 1]:length(d18Oc[, 1]), 1:Drange_pos[j, 2]), 2]) - min(d18Oc[c(Drange_pos[j, 1]:length(d18Oc[, 1]), 1:Drange_pos[j, 2]), 2])) / 2
                    t_SD[j] <- length(c(Drange_pos[j, 1]:length(d18Oc[, 1]), 1:Drange_pos[j, 2])) / 2
                    gr_SD[j] <- (max(GR[c(Drange_pos[j, 1]:length(GR[, 1]), 1:Drange_pos[j, 2]), 2]) - min(GR[c(Drange_pos[j, 1]:length(GR[, 1]), 1:Drange_pos[j, 2]), 2])) / 2
                    Tmod_SD[j] <- (max(SST[c(Drange_pos[j, 1]:length(SST[, 1]), 1:Drange_pos[j, 2]), 2]) - min(SST[c(Drange_pos[j, 1]:length(SST[, 1]), 1:Drange_pos[j, 2]), 2])) / 2
                }else{ # If the D range is contained within one year, range of d18Oc values is located between the lower and upper boundaries of D (+/- 1 SD)
                    Omod_SD[j] <- (max(d18Oc[Drange_pos[j, 1]:Drange_pos[j, 2], 2]) - min(d18Oc[Drange_pos[j, 1]:Drange_pos[j, 2], 2])) / 2
                    t_SD[j] <- length(Drange_pos[j, 1]:Drange_pos[j, 2]) / 2
                    gr_SD[j] <- (max(GR[Drange_pos[j, 1]:Drange_pos[j, 2], 2]) - min(GR[Drange_pos[j, 1]:Drange_pos[j, 2], 2])) / 2
                    Tmod_SD[j] <- (max(SST[Drange_pos[j, 1]:Drange_pos[j, 2], 2]) - min(SST[Drange_pos[j, 1]:Drange_pos[j, 2], 2])) / 2
                }
            }
            resmat <- cbind(resmat, Omod_SD, t_SD, gr_SD, Tmod_SD)
        }else{
            resmat <- cbind(resmat, matrix(rep(NA, 4 * length(resmat[,1])), ncol = 4)) # Add NA-filled columns if no MC error propagation is done
        }
        colnames(resmat) <- c("Dsam",
            "Osam",
            "Omod",
            "Residuals",
            "Time_of_year",
            "Instantaneous_growth_rate",
            "Modelled_temperature",
            "Omod_SD",
            "Time_of_Year_SD",
            "Instantaneous_growth_rate_SD",
            "Modelled_temperature_SD")

        if(plot == TRUE){
            plotdf <- as.data.frame(cbind(resmat, D_err, O_err))
            dev.new()
            ggplot2::ggplot(plotdf, ggplot2::aes(x = Dsam, y = Osam)) +
                ggplot2::geom_point() + # Plot original data
                ggplot2::geom_errorbar(ggplot2::aes(ymin = Osam - O_err, # Add error bars on measurement (1 SD)
                    ymax = Osam + O_err),
                    width = 100,
                    col = "black") +
                ggplot2::geom_errorbarh(ggplot2::aes(xmin = Dsam - D_err,
                    xmax = Dsam + D_err),
                    height = 0.05,
                    col = "black") +
                ggplot2::geom_point(ggplot2::aes(Dsam, Omod), col="red") + # Plot modelled d18Oc on top of data
                ggplot2::geom_errorbar(ggplot2::aes(ymin = Omod - Omod_SD, # Add error bars on model result (1 SD)
                    ymax = Omod + Omod_SD),
                    width = 100,
                    col = "red") +
                ggplot2::geom_errorbarh(ggplot2::aes(Dsam,
                    Omod,
                    xmin = Dsam - D_err,
                    xmax = Dsam + D_err),
                    height = 0.05,
                    col = "red")
        }
        return(resmat)
    }
}