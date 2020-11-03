#' Full ShellChron workflow wrapped in a single function
#' 
#' Takes starting parameters and names of input files
#' and directory and runs through all the steps of the
#' ShellChron model. Function includes options for plotting
#' and exporting raw data, which are parsed into underlying
#' formulae.
#' @param path String containing the path to the directory
#' that contains the input data.
#' @param file_name Name of the file that contains d18O data
#' @param mineral String containing the name of the 
#' mineralogy (either \code{"calcite"} or \code{"aragonite"})
#' Defaults to calcite.
#' @param t_int Time interval (in days; default = 1)
#' @param T_per Period of SST sinusoid (in days;
#' default = 365)
#' @param d18Ow Either a single value (constant d18Ow)
#'  or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about
#' seasonality in d18Ow. Defaults to constant d18Ow of
#' 0 permille VSMOW (the modern mean ocean value)
#' @param t_maxtemp Timing of the warmest day of the year
#' (in julian day; default = 182.5, or May 26th halfway
#' through the year)
#' @param MC Number of Monte Carlo simulations to apply for
#' error propagation. Default = 1000
#' @param plot Should an overview of the results of modelling
#' be plotted? \code{TRUE/FALSE}
#' @param plot_export Should the overview plot be exported as
#' a PDF file? \code{TRUE/FALSE}
#' @param export_raw Export tables containing all raw model
#' results before being merged into tidy tables? \code{TRUE/FALSE}
#' @param export_path Path where result files are exported
#' @return CSV tables of model results in the current working
#' directory, optional plots in PDF format and list object of
#' model results for further processing in the R workspace.
#' @references function dependencies: data_import, run_model, cumulative_day,
#' export_results
#' @importFrom grDevices dev.new dev.off pdf
#' @importFrom graphics lines points
#' @importFrom stats D approx lm loess pf qt rnorm sd spectrum ts weighted.mean
#'     window
#' @importFrom utils capture.output head read.csv tail write.csv
#' @examples
#' # find attached dummy data
#' example <- wrap_function(path = getwd(),
#'     file_name = system.file("extdata", "Virtual_shell.csv",
#'     package = "ShellChron"),
#'     "calcite",
#'     1,
#'     365,
#'     d18Ow = 0,
#'     t_maxtemp = 182.5,
#'     MC = 1000,
#'     plot = FALSE,
#'     plot_export = FALSE,
#'     export_raw = FALSE,
#'     export_path = tempdir()) # Run function
#' @export
wrap_function <- function(path, # Wrapping function for the entire model package
    file_name, # Give file name (don't forget to add the extention, should be in CSV format)
    mineral = "calcite", # Set mineralogy of the record, default is calcite. Aragonite is also supported
    t_int = 1, # Set time interval in days
    T_per = 365, # Set annual time period in days (default = 365)
    d18Ow = "default", # Set d18Ow value or vector (default = constant year-round at 0 VSMOW). Alternative options are either one value (assumed constant year-round) or a vector with length T_per / t_int and interval t_int specifying d18Ow evolution through one year.
    t_maxtemp = 182.5, # Define the day of the year at which temperature is heighest. Default = Assume that the day of maximum temperature is helfway through the year
    MC = 1000, # Number of MC simulations to include measurement error into error analysis. Default = 1000 (if MC = 0, error on D and d18O measurements not considered)
    plot = TRUE, # Should intermediate plots be given to track progress? WARNING: plotting makes the script much slower, especially for long datasets.
    plot_export = TRUE, # Should a plot of the results be saved as PDF?
    export_raw = FALSE, # Should the results of all individual model runs be exported as CSV files?
    export_path # Path where result files are exported
    ){

    # STEP 1: Import data
    setwd(path)
    importlist <- data_import(file_name)
    dat <- importlist[[1]]
    dynwindow <- importlist[[2]]
    G_per <- T_per # Period of growth rate sinusoid should equal that of the temperature sinusoid (which is given)

    # STEP 2: Run the model
    resultlist <- run_model(dat, dynwindow, mineral, d18Ow, T_per, G_per, t_int, t_maxtemp, MC, plot = TRUE)
    resultarray <- resultlist[[1]]
    parmat <- resultlist[[2]]
    
    # STEP 3: Align model results to cumulative timescale
    print("Calculating cumulative day of the year results...")
    suppressWarnings(resultarray[, , 3] <- cumulative_day(resultarray, TRUE, TRUE, export_path)) # Calculate cumulative day of the year for all model runs and replace matrix in result array
    
    # STEP 4: Order and export results and statistics
    export_results(export_path, dat, resultarray, parmat, MC, dynwindow, plot, plot_export, export_raw) # Export results of model
    return(resultlist)
}