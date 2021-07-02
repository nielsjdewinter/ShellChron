#' Function that carries out a sinusoidal regression
#' 
#' Fits a sinusoid through data provided as an \code{x} and \code{y} 
#' vector and returns a list containing both the fitted curve and the
#' parameters of that curve.
#' Used to produce initial values for modeling data windows and later
#' to find peaks in modeled julian day values to align the result to
#' a cumulative age timeline.
#' 
#' @param x Vector of \code{x} values of input data
#' @param y Vector of \code{y} values of input data
#' @param fixed_period Optional variable for fixing the period of the sinusoid 
#' in the depth domain. Defaults to \code{NA}, period is not fixed. Supply a 
#' single value to fix the period.
#' @param plot Should the fitting result be plotted? \code{TRUE/FALSE}
#' @return A list containing a vector of parameters of the fitted sinusoid
#' and the fitted values belonging to each \code{x} value.
#' Fitting parameters:
#' \code{I} = the mean annual value of the sinusoid (height)
#' \code{A} = the amplitude of the sinusoid
#' \code{Dper} = the period of the sinusoid in \code{x} domain
#' \code{peak} = the location of the peak in the sinusoid
#' \code{R2adj} = the adjusted \code{R^2} value of the fit
#' \code{p} = the p-value of the fit
#' @examples
#' # Create dummy data
#' x <- seq(1000, 11000, 1000)
#' y <- sin((2 * pi * (seq(1, 11, 1) - 8 + 7 / 4)) / 7)
#' sinlist <- sinreg(x, y, plot = FALSE) # Run the function
#' @export
sinreg <- function(x, # Function to perform sinusoid regression meant to estimate the starting parameters for fitting growth and temperature sinusoids
    y,
    fixed_period = NA, # Option to fix the period in depth domain. Default = NA (Period not fixed)
    plot = FALSE # Plot results?
    ){

    if(is.na(fixed_period)){ # If period is not fixed, run simple periodogram to find the most likely period
        Ots <- ts(y) # Turn y into a time series
        ssp <- spectrum(Ots, plot = FALSE) # Create periodogram of y
        Nper <- 1/ssp$freq[ssp$spec==max(ssp$spec)] # Estimate period in terms of sample number
        Dper <- Nper * diff(range(x))/length(x) # Convert period to depth domain
    }else if(is.numeric(fixed_period) & length(fixed_period) == 1){ # Check if a valid entry is provided for Dper
        Dper <- fixed_period
    }else{
        stop("ERROR: Supplied value for fixed period is not recognized")
    }

    sinlm <- lm(y ~ sin(2*pi/Dper*x)+cos(2*pi/Dper*x)) # Run linear model, cutting up d18O = Asin(2*pi*D) into d18O = asin(D)+bcos(D)
    sinm<-sinlm$fitted # Extract model result
    if(plot == TRUE){
        dev.new()
        plot(x, y); lines(x, sinm, col = "red")
    }

    coeff<-summary(sinlm)[["coefficients"]] # Extract a, b and intercept + uncertainties

    # Calculate coefficients of the form d18O = I + Asin(2*pi*D + p)
    I<-coeff[1,1] # Intercept (mean annual value)
    A<-sqrt(coeff[2,1]^2+coeff[3,1]^2) # Amplitude of seasonality
    phase<--acos(coeff[2,1]/A) # Phase of sinusoid
    peak<-(0.25+(phase/(2*pi))) * Dper # timing of seasonal peak
    R2adj<-summary(sinlm)$adj.r.squared # Goodness of model fit (adjusted R2)
    p <- as.numeric(pf(summary(sinlm)$fstatistic[1],summary(sinlm)$fstatistic[2],summary(sinlm)$fstatistic[3],lower.tail=F)) # p-value of model

    return(list(c(I,A,Dper,peak,R2adj,p),sinm)) # Return results of regression
}