#' Function that propagates measurement uncertainty through model results
#' 
#' Function to propagate combined errors on \code{x} (= \code{Dsam}) and
#' \code{y} (= \code{Osam}) on the modelled X (= \code{D}) and Y 
#' \code{d18Oc} values by means of direct projection of y-uncertainty
#' on \code{x} and then combine the errors on both in the \code{x} domain
#' 
#' Note: projection y_err on x_err leads to large X errors on shallow
#' slopes due to numerical calculation of fist derivative.
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
#' Y <- cbind(seq(1, 39, 1), 0.9 * sin((2 * pi * (seq(1, 39, 1) - 9 +
#'     25 / 4)) / 25))
#' # Run function
#' result <- mc_err_proj(x, x_err, y, y_err, X, Y, 1000)
#' @export
mc_err_proj <- function(x, x_err, y, y_err, X, Y, MC){ # Function to propagate
# combined errors on x and y on the modelled X and Y values by means of local
# projection of y uncertainty on x and subsequent combination of uncertainties
# in X domain
    dYdX <- diff(Y[, 2]) / diff(X) # Create first derivative of Y by X
    dYdX[which(abs(dYdX) <= 10 ^ (floor(log(mean(abs(dYdX)), 10)) - 1))] <-
        sign(dYdX[which(abs(dYdX) <= 10 ^ (floor(log(mean(abs(dYdX)),
            10)) - 1))]) * 10 ^ (floor(log(mean(abs(dYdX)), 10)) - 1) # Remove
            # small absolute values (more than one order of magnitude smaller
            # than the mean), preserving their sign
    dYdX <- append(dYdX, dYdX[length(dYdX)]) # repeat last value to increase the
    # length of the vector to match X/Y (366 days)
    localslope <- dYdX[apply(abs(outer(x, X, FUN = "-")), 1, which.min)] # Find
    # the local slope belonging to each sample position
    x_err_proj <- y_err / localslope # Project uncertainty on Y on x domain
    # using slope
    x_err_comb <- sqrt(x_err_proj ^2 + x_err ^2) # Combine uncertainties on x
    # and y in the X domain
    return(x_err_comb)
}