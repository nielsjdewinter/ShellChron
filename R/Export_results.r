#' Function to merge and export the results of the ShellChron model
#' 
#' Takes the input data and model results and reformats
#' them to tables of key parameters such as growth rate
#' and shell age for each datapoint for easy plotting.
#' This final function also combines uncertainties in the
#' model result arising from uncertainties in input data
#' (provided by the user) and uncertainties of the model
#' (from overlapping modelling windows).
#' Includes some optional plotting options.
#' @param path Path where result files are exported
#' @param dat Matrix containing the input data
#' @param resultarray Array containing the full results of
#' the optimized growth model
#' @param parmat Matrix listing all optimized growth
#' rate and SST parameters used to model d18O in each data
#' window
#' @param MC Number of Monte Carlo simulations to apply for
#' error propagation. Default = 1000
#' @param dynwindow Information on the position and length 
#' of modelling windows
#' @param plot Should an overview of the results of modelling
#' be plotted? \code{TRUE/FALSE}
#' @param plot_export Should the overview plot be exported as
#' a PDF file? \code{TRUE/FALSE}
#' @param export_raw Export tables containing all raw model
#' results before being merged into tidy tables? \code{TRUE/FALSE}
#' @return CSV tables of model results in the current working
#' directory + optional plots in PDF format
#' @references package dependencies: tidyverse 1.3.0; ggpubr 0.4.0; magrittr  
#' function dependencies: sd_wt
#' @examples
#' # Create dummy input data column by column
#' dat <- as.data.frame(seq(1000, 40000, 1000))
#' colnames(dat) <- "D"
#' dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
#' dat$YEARMARKER <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
#' dat$D_err <- rep(100, 40)
#' dat$d18Oc_err <- rep(0.1, 40)
#' 
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
#' # Create dummy seasonality data
#' seas <- as.data.frame(seq(1, 365, 1))
#' colnames(seas) <- "t"
#' seas$SST <- 15 + 10 * sin((2 * pi * (seq(1, 365, 1) - 182.5 +
#'     365 / 4)) / 365)
#' seas$GR <- 10 + 10 * sin((2 * pi * (seq(1, 365, 1) - 100 + 365 / 4)) / 365)
#' seas$d18O <- (exp((18.03 * 1000 / (seas$SST + 273.15) - 32.42) / 1000) - 1) *
#'     1000 + (0.97002 * 0 - 29.98)
#' # Apply dummy seasonality data to generate other tabs of testarray
#' testarray[, , 1] <- seas$d18O[match(testarray[, , 3], seas$t)] # d18O values
#' tab <- testarray[, , 1]
#' tab[which(!is.na(tab))] <- 0.1
#' testarray[, , 2] <- tab # dummy d18O residuals
#' testarray[, , 4] <- seas$GR[match(testarray[, , 3], seas$t)] # growth rates
#' testarray[, , 5] <- seas$SST[match(testarray[, , 3], seas$t)] # temperature
#' tab[which(!is.na(tab))] <- 0.1
#' testarray[, , 6] <- tab # dummy d18O SD
#' tab[which(!is.na(tab))] <- 20
#' testarray[, , 7] <- tab # dummy time SD
#' tab[which(!is.na(tab))] <- 3
#' testarray[, , 8] <- tab # dummy GR SD
#' tab[which(!is.na(tab))] <- 1
#' testarray[, , 9] <- tab # dummy temperature SD
#' darray <- array(rep(as.matrix(dat), 9), dim = c(40, 5, 9))
#' testarray[, 1:5, ] <- darray
#' 
#' # Create dummy dynwindow data
#' dynwindow <- as.data.frame(seq(1, 31, 1))
#' colnames(dynwindow) <- "x"
#' dynwindow$y <- rep(10, 31)
#' 
#' dimnames(testarray) <- list(
#'     paste("sample", 1:length(testarray[, 1, 3])),
#'     c(colnames(dat), paste("window", 1:length(dynwindow$x))),
#'     c("Modelled_d18O",
#'         "d18O_residuals",
#'         "Time_of_year",
#'         "Instantaneous_growth_rate",
#'         "Modelled temperature",
#'         "Modelled_d18O_SD",
#'         "Time_of_Year_SD",
#'         "Instantaneous_growth_rate_SD",
#'         "Modelled_temperature_SD")
#' )
#' 
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
#' parsSD <- c(3, 10, 3, 5, 10, 3, 5) # Artificial variability in parameters
#' parmat <- matrix(rnorm(length(pars) * length(dynwindow$x)), nrow =
#'     length(pars)) * parsSD + matrix(rep(pars, length(dynwindow$x)),
#'     nrow = length(pars))
#' rownames(parmat) <- c("T_amp", "T_pha", "T_av", "G_amp", "G_pha", "G_av",
#'     "G_skw")
#' # Run export function
#' test <- export_results(path = tempdir(),
#'     dat,
#'     testarray,
#'     parmat,
#'     MC = 1000,
#'     dynwindow,
#'     plot = FALSE,
#'     plot_export = FALSE,
#'     export_raw = FALSE)
#' @export
export_results <- function(path, # Path where result files are exported
    dat, # raw data
    resultarray, # Array containing all model results 
    parmat, # matrix of parameters per window
    MC, # Include number of simulations just for error verification (if MC > 0, errors are included in the export)
    dynwindow, # Include the size of the windows for pooling standard deviations AND FOR ADDING WEIGHINGS TO STATISTICS BASED ON PLACE IN WINDOW
    plot = FALSE, # Create a result plot?
    plot_export = TRUE, # Export a result plot?
    export_raw = FALSE # Export all the raw results of the model (of individual windows)?
    ){
    
    Day <- sd.day <- N <- se.day <- d18O_mod <- sd.d18O_mod <- se.d18O_mod <-
        GR <- sd.GR <- se.GR <- SST <- sd.SST <- se.SST <- parameter <- 
        par_value <- stdev <- se.pars <- SD <- d18Oc <- mean.day <- CL95.day <-
        mean.d18O_mod <- CL95.d18O_mod <- mean.GR <- CL95.GR <- NULL # Predefine variables to circumvent global variable binding error
    
    # Define weights to give more priority to datapoints in the center of the modelling window than those on the edge
    weights <- matrix(NA, ncol = length(dynwindow$x), nrow = length(dynwindow$x) + dynwindow$y[length(dynwindow$x)] - 1) # Create template matrix
    for(i in 1:length(dynwindow$x)){ # Loop through matrix and add weights for each position in the resultarray that contains a value
        weights[dynwindow$x[[i]]:(dynwindow$x[[i]] + dynwindow$y[[i]] - 1), i] <- dynwindow$y[[i]] / 2 - abs(dynwindow$x[[i]]:(dynwindow$x[[i]] + dynwindow$y[[i]] - 1) - (dynwindow$x[[i]] + (dynwindow$y[[i]] - 1) / 2))
    }
    weights <- cbind(resultarray[, 1:5, 3], weights)
    weightstidy <- tidyr::gather(as.data.frame(weights), "window", "weight", (length(dat[1, ]) + 1):ncol(weights), factor_key = TRUE) # Convert weights to Tidy data for plotting

    JDtidy <- tidyr::gather(as.data.frame(resultarray[, , 3]), "window", "Day", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE) # Convert modelled time results to Tidy data for plotting
    JDtidy$weights <- weightstidy$weight # Add weights to JDtidy

    JDstats <- JDtidy %>% # Summarize modelled time statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            mean.day = mean(Day, na.rm = TRUE),  # Calculate means per sample
            sd.day = sd_wt(Day, weights, na.rm = TRUE),  # Calculate stdevs per sample
            N = dplyr::n_distinct(Day, na.rm = TRUE), # Calculate the number of modelled values, excluding NA's
            se.day = sd.day / sqrt(N), # Calculate the standard error
            CL95.day = qt(0.95, N) * se.day # Calculate the 95% confidence level
        )
    JDstats$sd.day[which(JDstats$N == 1)] <- NaN

    d18Otidy <- tidyr::gather(as.data.frame(resultarray[, , 1]), "window", "d18O_mod", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE) # Convert modelled d18O results to Tidy data for plotting
    d18Otidy$weights <- weightstidy$weight # Add weights to d18Otidy

    d18Ostats <- d18Otidy %>% # Summarize modelled d18O statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            mean.d18O_mod = mean(d18O_mod, na.rm = TRUE),  # Calculate means per sample
            sd.d18O_mod = sd_wt(d18O_mod, weights, na.rm = TRUE),  # Calculate stdevs per sample
            N = dplyr::n_distinct(d18O_mod, na.rm = TRUE), # Calculate the number of modelled values, excluding NA's
            se.d18O_mod = sd.d18O_mod / sqrt(N), # Calculate the standard error
            CL95.d18O_mod = qt(0.95, N) * se.d18O_mod # Calculate the 95% confidence level
        )
    d18Ostats$sd.d18O_mod[which(d18Ostats$N == 1)] <- NaN

    GRtidy <- tidyr::gather(as.data.frame(resultarray[, , 4]), "window", "GR", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE) # Convert modelled growth rate results to Tidy data for plotting
    GRtidy$weights <- weightstidy$weight # Add weights to GRtidy

    GRstats <- GRtidy %>% # Summarize modelled growth rate statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            mean.GR = mean(GR[GR>0.1], na.rm = TRUE),  # Calculate means per sample, excluding NA's and instances where growth rate is near-zero
            sd.GR = sd_wt(GR[GR>0.1], weights, na.rm = TRUE),  # Calculate stdevs per sample, excluding NA's and instances where growth rate is near-zero
            N = dplyr::n_distinct(GR[GR>0.1], na.rm = TRUE), # Calculate the number of modelled values, excluding NA's and instances where growth rate is near-zero
            se.GR = sd.GR / sqrt(N), # Calculate the standard error
            CL95.GR = qt(0.95, N) * se.GR # Calculate the 95% confidence level
        )
    GRstats$sd.GR[which(GRstats$N == 1)] <- NaN

    Ttidy <- tidyr::gather(as.data.frame(resultarray[, , 5]), "window", "SST", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE) # Convert modelled temperature results to Tidy data for plotting
    Ttidy$weights <- weightstidy$weight # Add weights to Ttidy

    Tstats <- Ttidy %>% # Summarize modelled growth rate statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            mean.SST = weighted.mean(SST[SST>0.1], na.rm = TRUE),  # Calculate means per sample, excluding NA's and instances where SSTowth rate is near-zero
            sd.SST = sd_wt(SST[SST>0.1], weights, na.rm = TRUE),  # Calculate stdevs per sample, excluding NA's and instances where SSTowth rate is near-zero
            N = dplyr::n_distinct(SST[SST>0.1], na.rm = TRUE), # Calculate the number of modelled values, excluding NA's and instances where SSTowth rate is near-zero
            se.SST = sd.SST / sqrt(N), # Calculate the standard error
            CL95.SST = qt(0.95, N) * se.SST # Calculate the 95% confidence level
        )
    Tstats$sd.SST[which(Tstats$N == 1)] <- NaN

    parmat2 <- data.frame(rownames(parmat), parmat)
    colnames(parmat2)[1] <- "parameter"
    partidy <- tidyr::gather(parmat2, "window", "par_value", 2:length(parmat2[1,]), factor_key = TRUE)

    parstats <- partidy %>% # Summarize model parameters
        ggpubr::group_by(parameter) %>%
        dplyr::summarize(
            means = mean(par_value), # Calculate means per parameter
            stdev = sd(par_value), # Calculate standard deviation per parameter
            N = dplyr::n(), # Count number of modelled values per parameter (= equal to number of windows)
            se.pars = stdev / sqrt(N), # Calculate standard error
            CL95 = qt(0.95, N) * se.pars
        )
    
    if(MC > 0){
        print("Recalculating export statistics by including propagated uncertainties")
        # Include errors propagated from those on D and d18Oc data into the statistics

        # Propagate errors on modelled d18O
        d18Otidy_err <- d18Otidy
        d18Otidy_err$SD <- tidyr::gather(as.data.frame(resultarray[, , 6]), "window", "d18O", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE)$d18O # Convert modelled time errors to Tidy data for plotting
        d18Otidy_err$N <- dynwindow$y[as.numeric(d18Otidy$window)] # Add window size for calculating pooled SD
        d18Otidy_err <- d18Otidy_err[-which(is.na(d18Otidy_err$d18O_mod)), ] # Remove empty cells in matrix
        d18Otidy_err$SD[which(d18Otidy_err$SD == 0)] <- min(d18Otidy_err$SD[which(d18Otidy_err$SD > 0)]) # Replace zeroes with smallest SD to prevent division by zero

        d18Ostats2 <- d18Otidy_err %>% # Summarize modelled d18O statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            weighted.mean.d18O = weighted.mean(d18O_mod, 1 / SD ^ 2 * weights, na.rm = TRUE),  # Calculate weighted means per sample
            pooled.sd.d18O = sqrt(sum(SD ^ 2 * (N - 1) * weights, na.rm = TRUE) / ((sum(N, na.rm = TRUE) - dplyr::n()) * mean(weights))) # Calculate pooled standard deviation resulting from error propagations and the weighted mean of the variances taking weights derived from position in the window into account
        )

        # Aggregate propagated errors into statistics matrices
        d18Ostats$mean.d18O_mod <- d18Ostats2$weighted.mean.d18O # Replace means by weighed means, taking into account the propagated error on individual estimates
        d18Ostats$sd.d18O_mod <- sqrt(d18Ostats$sd.d18O_mod ^ 2 + d18Ostats2$pooled.sd.d18O ^ 2) # Combine errors from the model and the errors on input
        d18Ostats$se.d18O_mod <- d18Ostats$sd.d18O_mod / sqrt(d18Ostats$N) # Propagate new errors onto standard error
        d18Ostats$CL95.d18O_mod <- qt(0.95, d18Ostats$N) * d18Ostats$se.d18O_mod # Propagate new errors onto confidence interval


        # Propagate errors on Time of Day calculations
        JDtidy_err <- JDtidy
        JDtidy_err$SD <- tidyr::gather(as.data.frame(resultarray[, , 7]), "window", "Day", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE)$Day # Convert modelled time errors to Tidy data for plotting
        JDtidy_err$N <- dynwindow$y[as.numeric(JDtidy$window)] # Add window size for calculating pooled SD
        JDtidy_err <- JDtidy_err[-which(is.na(JDtidy_err$Day)), ] # Remove empty cells in matrix
        JDtidy_err$SD[which(JDtidy_err$SD == 0)] <- min(JDtidy_err$SD[which(JDtidy_err$SD > 0)]) # Replace zeroes with smallest SD to prevent division by zero

        JDstats2 <- JDtidy_err %>% # Summarize modelled JD statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            weighted.mean.day = weighted.mean(Day, 1 / SD ^ 2 * weights, na.rm = TRUE),  # Calculate weighted means per sample
            pooled.sd.day = sqrt(sum(SD ^ 2 * (N - 1) * weights, na.rm = TRUE) / ((sum(N, na.rm = TRUE) - dplyr::n()) * mean(weights))) # Calculate pooled standard deviation resulting from error propagations and the weighted mean of the variances taking weights derived from position in the window into account
        )

        # Aggregate propagated errors into statistics matrices
        JDstats$mean.day <- JDstats2$weighted.mean.day # Replace means by weighed means, taking into account the propagated error on individual estimates
        JDstats$sd.day <- sqrt(JDstats$sd.day ^ 2 + JDstats2$pooled.sd.day ^2) # Combine errors from the model and the errors on input
        JDstats$se.day <- JDstats$sd.day / sqrt(JDstats$N) # Propagate new errors onto standard error
        JDstats$CL95.day <- qt(0.95, JDstats$N) * JDstats$se.day # Propagate new errors onto confidence interval


        # Propagate errors on modelled growth rate
        GRtidy_err <- GRtidy
        GRtidy_err$SD <- tidyr::gather(as.data.frame(resultarray[, , 8]), "window", "GR", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE)$GR # Convert modelled time errors to Tidy data for plotting
        GRtidy_err$N <- dynwindow$y[as.numeric(GRtidy$window)] # Add window size for calculating pooled SD
        GRtidy_err <- GRtidy_err[-which(is.na(GRtidy_err$GR)), ] # Remove empty cells in matrix
        GRtidy_err$SD[which(GRtidy_err$SD == 0)] <- min(GRtidy_err$SD[which(GRtidy_err$SD > 0)]) # Replace zeroes with smallest SD to prevent division by zero

        GRstats2 <- GRtidy_err %>% # Summarize modelled GR statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            weighted.mean.GR = weighted.mean(GR, 1 / SD ^ 2 * weights, na.rm = TRUE),  # Calculate weighted means per sample
            pooled.sd.GR = sqrt(sum(SD ^ 2 * (N - 1) * weights, na.rm = TRUE) / ((sum(N, na.rm = TRUE) - dplyr::n()) * mean(weights))) # Calculate pooled standard deviation resulting from error propagations and the weighted mean of the variances taking weights derived from position in the window into account
        )

        # Aggregate propagated errors into statistics matrices
        GRstats$mean.GR <- GRstats2$weighted.mean.GR # Replace means by weighed means, taking into account the propagated error on individual estimates
        GRstats$sd.GR <- sqrt(GRstats$sd.GR ^ 2 + GRstats2$pooled.sd.GR ^2) # Combine errors from the model and the errors on input
        GRstats$se.GR <- GRstats$sd.GR / sqrt(GRstats$N) # Propagate new errors onto standard error
        GRstats$CL95.GR <- qt(0.95, GRstats$N) * GRstats$se.GR # Propagate new errors onto confidence interval


        # Propagate errors on modelled temperature
        Ttidy_err <- Ttidy
        Ttidy_err$SD <- tidyr::gather(as.data.frame(resultarray[, , 9]), "window", "T", (length(dat[1, ]) + 1):length(resultarray[1, , 1]), factor_key = TRUE)$T # Convert modelled time errors to Tidy data for plotting
        Ttidy_err$N <- dynwindow$y[as.numeric(Ttidy$window)] # Add window size for calculating pooled SD
        Ttidy_err <- Ttidy_err[-which(is.na(Ttidy_err$SST)), ] # Remove empty cells in matrix
        Ttidy_err$SD[which(Ttidy_err$SD == 0)] <- min(Ttidy_err$SD[which(Ttidy_err$SD > 0)]) # Replace zeroes with smallest SD to prevent division by zero

        Tstats2 <- Ttidy_err %>% # Summarize modelled T statistics
        ggpubr::group_by(D) %>%
        dplyr::summarize(
            weighted.mean.SST = weighted.mean(SST, 1 / SD ^ 2 * weights, na.rm = TRUE),  # Calculate weighted means per sample
            pooled.sd.SST = sqrt(sum(SD ^ 2 * (N - 1) * weights, na.rm = TRUE) / ((sum(N, na.rm = TRUE) - dplyr::n()) * mean(weights))) # Calculate pooled standard deviation resulting from error propagations and the weighted mean of the variances taking weights derived from position in the window into account
        )

        # AgTegate propagated errors into statistics matrices
        Tstats$mean.SST <- Tstats2$weighted.mean.SST # Replace means by weighed means, taking into account the propagated error on individual estimates
        Tstats$sd.SST <- sqrt(Tstats$sd.SST ^ 2 + Tstats2$pooled.sd.SST ^2) # Combine errors from the model and the errors on input
        Tstats$se.SST <- Tstats$sd.SST / sqrt(Tstats$N) # Propagate new errors onto standard error
        Tstats$CL95.SST <- qt(0.95, Tstats$N) * Tstats$se.SST # Propagate new errors onto confidence interval
        print("Preparing plots")
    }

    if(plot == TRUE | plot_export == TRUE){ # Check if plots are needed
        # Create depth-time plot
        Dtplot <- ggplot2::ggplot(JDtidy, ggplot2::aes(D, Day)) +
            ggplot2::geom_point(ggplot2::aes(colour = d18Oc)) +
            ggplot2::scale_colour_gradient2(midpoint = mean(JDtidy$d18Oc)) +
            ggplot2::geom_line(data = JDstats, ggplot2::aes(D, mean.day), size = 1) +
            ggplot2::geom_line(data = JDstats, ggplot2::aes(D, mean.day + CL95.day), size = 1, alpha = 0.5) +
            ggplot2::geom_line(data = JDstats, ggplot2::aes(D, mean.day - CL95.day), size = 1, alpha = 0.5) +
            ggplot2::ggtitle("Plot of Height vs. Time") +
            ggplot2::xlab("Record length") +
            ggplot2::scale_y_continuous("Age (days)", seq(0, 365 * ceiling(max(JDstats$mean.day + JDstats$CL95.day, na.rm = TRUE) / 365), 365))

        # Create d18O plot
        d18Oplot <- ggplot2::ggplot(d18Otidy, ggplot2::aes(D, d18Oc)) +
            ggplot2::geom_point() +
            ggplot2::geom_line(data = d18Ostats, ggplot2::aes(D, mean.d18O_mod), size = 1) +
            ggplot2::geom_line(data = d18Ostats, ggplot2::aes(D, mean.d18O_mod + CL95.d18O_mod, alpha = 0.5, col = "darkblue"), size = 1) +
            ggplot2::geom_line(data = d18Ostats, ggplot2::aes(D, mean.d18O_mod - CL95.d18O_mod, alpha = 0.5, col = "darkred"), size = 1) +
            ggplot2::ggtitle("Plot of measured and modelled d18O vs. Record Length") +
            ggplot2::xlab("Record length") +
            ggplot2::ylab("d18O_carbonate") +
            ggplot2::theme(legend.position = "none") # Remove legend

        # Create growth rate plot

        GRplot <- ggplot2::ggplot(GRtidy, ggplot2::aes(D, GR)) +
            ggplot2::geom_point(ggplot2::aes(colour = d18Oc)) +
            ggplot2::scale_colour_gradient2(midpoint = mean(JDtidy$d18Oc)) +
            ggplot2::geom_line(data = GRstats, ggplot2::aes(D, mean.GR), size = 1) +
            ggplot2::geom_line(data = GRstats, ggplot2::aes(D, mean.GR + CL95.GR, alpha = 0.5), size = 1) +
            ggplot2::geom_line(data = GRstats, ggplot2::aes(D, mean.GR - CL95.GR, alpha = 0.5), size = 1) +
            ggplot2::ggtitle("Plot of modelled growth rate vs Record Length") +
            ggplot2::xlab("Record length") +
            ggplot2::ylab("Growth rate") +
            ggplot2::theme(legend.position = "none") # Remove legend

        Combined_plots <- ggpubr::ggarrange(Dtplot, d18Oplot, GRplot, labels = c("A", "B", "C"), ncol = 3, nrow = 1) # Combine plots

        if(plot == TRUE){
            dev.new()
            print(Combined_plots)
        }

        if(plot_export == TRUE){
            pdf("Model result plots.pdf", width = 30, height = 10)
            print(Combined_plots)
            dev.off()
        }
    }
    print("Start exporting files to directory")
    if(export_raw == TRUE){
    # Write away all raw results of modelling
        d18Oraw_p <- file.path(path, "modelled_d18O_raw.csv")
        write.csv(resultarray[, , 1], d18Oraw_p)
        resraw_p <- file.path(path, "residuals_raw.csv")
        write.csv(resultarray[, , 2], resraw_p)
        JDraw_p <- file.path(path, "Day_of_year_raw.csv")
        write.csv(resultarray[, , 3], JDraw_p)
        IGRraw_p <- file.path(path, "Instantaneous_growth_rate_raw.csv")
        write.csv(resultarray[, , 4], IGRraw_p)
        SSTraw_p <- file.path(path, "SST_raw.csv")
        write.csv(resultarray[, , 5], SSTraw_p)
        d18OSDraw_p <- file.path(path, "Modelled_d18O_SD_raw.csv")
        write.csv(resultarray[, , 6], d18OSDraw_p)
        JDSDraw_p <- file.path(path, "Day_of_Year_SD_raw.csv")
        write.csv(resultarray[, , 7], JDSDraw_p)
        IGRSDraw_p <- file.path(path, "Instantaneous_growth_rate_SD_raw.csv")
        write.csv(resultarray[, , 8], IGRSDraw_p)
        SSTSDraw_p <- file.path(path, "SST_SD_raw.csv")
        write.csv(resultarray[, , 9], SSTSDraw_p)
        parraw_p <- file.path(path, "modelled_parameters_raw.csv")
        write.csv(parmat, parraw_p)
    }

    # Write avay summary statistics of modelling
    AMR_p <- file.path(path, "Age_model_results.csv")
    write.csv(JDstats, AMR_p)
    d18OR_p <- file.path(path, "d18O_model_results.csv")
    write.csv(d18Ostats, d18OR_p)
    GRR_p <- file.path(path, "Growth_rate_results.csv")
    write.csv(GRstats, GRR_p)
    SSTR_p <- file.path(path, "SST_results.csv")
    write.csv(Tstats, SSTR_p)
    parR_p <- file.path(path, "Model_parameter_results.csv")
    write.csv(parstats, parR_p)
    print("DONE!")}