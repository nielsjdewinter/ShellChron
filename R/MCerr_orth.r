#' Function that propagates measurement uncertainty through model results
#' 
#' Function to propagate combined errors on \code{x} (= \code{Dsam}) and
#' \code{y} (= \code{Osam}) on the modelled X (= \code{D}) and Y 
#' \code{d18Oc} values by means of orthogonal projection of uncertainty
#' on \code{x} and \code{y} onto the model curve
#' @param x Vector of \code{x} values of input data
#' @param x_err Vector of uncertainties on \code{x} values
#' @param y Vector of \code{y} values of input data
#' @param y_err Vector of uncertainties on \code{y} values
#' @param X Vector of modelled \code{X} values on which the uncertainty is
#' to be projected
#' @param Y Matrix of modelled x and \code{Y} values
#' @param MC Number of Monte Carlo simulations to apply for error propagation
#' Default = 1000
#' @return A vector listing the standard deviations of propagated errors 
#' propagated on all \code{X} values.
#' @examples
#' # Create dummy data for input
#' x <-  seq(1, 40, 1)
#' x_err <- rep(0.1, 40)
#' y <- sin((2 * pi * (seq(1, 40, 1) - 8 + 30 / 4)) / 30)
#' y_err <- rep(0.1, 40)
#' X <- seq(1.5, 39.5, 1)
#' Y <- cbind(seq(1, 39, 1), 0.9 * sin((2 * pi * (seq(1, 39, 1) - 9 + 25 / 4)) / 25))
#' # Run function
#' result <- MCerr_orth(x, x_err, y, y_err, X, Y, 1000)
#' @export
MCerr_orth <- function(x, x_err, y, y_err, X, Y, MC){ # Function to propagate combined errors on x and y on the modelled X and Y values by means of orthogonal projection of x and y uncertainty on the model curve
    xmat <- matrix(rnorm(MC * length(x)), nrow = length(x)) * x_err + matrix(rep(x, MC), nrow = length(x)) # Create matrix of simulated X values
    ymat <- matrix(rnorm(MC * length(y)), nrow = length(y)) * y_err + matrix(rep(y, MC), nrow = length(y)) # Create matrix of simulated Y values
    Xarray <- sqrt( # create array of the length of vectors connecting each MC simulated x-y pair and each modelled X-Y pair.
        ((outer(xmat, X, FUN = "-") - mean(outer(xmat, X, FUN = "-"))) / sd(outer(xmat, X, FUN = "-"))) ^ 2 + # Square of the normalized difference between MC-simulated X values and modelled D
        ((outer(ymat, Y[, 2], FUN = "-") - mean(outer(ymat, Y[, 2], FUN = "-"))) / sd(outer(ymat, Y[, 2], FUN = "-"))) ^ 2 # Square of the normalized difference between MC-simulated y values and modelled Y values
    )
    posmat <- apply(Xarray, c(1,2), which.min) # For each MC-simulated x-y pair, find the position of the closest model value
    Xsimmat <- matrix(X[posmat], nrow = length(x)) # Find the X value that belongs to each position in posmat
    X_err_comb <- apply(Xsimmat, 1, sd) # calculate the standard deviation in X domain resulting from orthogonal projection of errors on X and Y from variability within the X values
    return(X_err_comb)
}