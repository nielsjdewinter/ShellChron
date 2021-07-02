#' Function that propagates measurement uncertainty through model results
#' 
#' Function to propagate combined errors on \code{x} (= \code{Dsam}) and
#' \code{y} (= \code{Osam}) on the modeled X (= \code{D}) and Y 
#' \code{d18Oc} values by means of projection of uncertainties
#' through the modeled \code{X-Y} relationship
#'
#' Note: projection leads to large uncertainties on shallow parts of the
#' \code{X-Y} curve
#' @param x Vector of \code{x} values of input data
#' @param x_err Vector of uncertainties on \code{x} values
#' @param y Vector of \code{y} values of input data
#' @param y_err Vector of uncertainties on \code{y} values
#' @param X Vector of modeled \code{X} values on which the uncertainty is
#' to be projected
#' @param Y Matrix of modeled x and \code{Y} values
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
#' result <- mc_err_form(x, x_err, y, y_err, X, Y, 1000)
#' @export
mc_err_form <- function(x,
    x_err,
    y,
    y_err,
    X,
    Y,
    MC = 1000){
    
    xmat <- matrix(rnorm(MC * length(x)), nrow = length(x)) * x_err + matrix(rep(x, MC), nrow = length(x)) # Create matrix of simulated X values
    Xpos <- apply(abs(outer(xmat, X, FUN = "-")), c(1,2), which.min) # find closest position in X vector (day) for each simulated X value
    Xpos_stat <- cbind(rowMeans(Xpos), apply(Xpos, 1, sd)) # Find mean and standard deviation of positions in X vector (day) for each sample
    
    ymat <- matrix(rnorm(MC * length(y)), nrow = length(y)) * y_err + matrix(rep(Y[Xpos_stat[, 1], 2], MC), nrow = length(y)) # Create matrix of Monte Carlo simulated Y values projected on the X-Y curve
    Xpos_mat <- outer(round(Xpos_stat[, 1]), seq(-20, 20, 1), "+")  %% length(D) + 1 # Create matrix of D positions around the mean position for each sample to match with simulated Y values. Window is +/- 20 positions
    Xpos_matO <- matrix(Y[Xpos_mat, 2], nrow = length(x)) # Find Y values for each position in Xpos_mat
    Xpos_matO2 <- apply(outer(Xpos_matO, ymat, FUN = "-"), c(2,4), diag) # Match and subtract each simulated Y value with the local environment of Y values in the X-Y curve
    Ypos <- apply(abs(Xpos_matO2), c(1, 3), which.min) + Xpos_stat[, 1] # Find the position of the closest Y value in the Y window (Xpos_matO) for each Y simulation
    Ypos_stat <- cbind(rowMeans(Ypos), apply(Ypos, 1, sd)) # Find mean and standard deviation of positions in X vector (day) for each sample
    Ypos_stat[which(Ypos_stat[, 2] == 0), 2] <- sd(seq(-20, 20, 1)) # If SD of Y falls outside the window of +/- 20 positions SD is assumed to be equal to that of the uniform distribution
    
    pos_err_comb <- sqrt(Ypos_stat[, 2] ^ 2 + Xpos_stat[, 2] ^ 2) # Find combined error on position (day)
    Xpos_minmax <- cbind(Xpos_stat[, 1] - pos_err_comb, Xpos_stat[, 1] + pos_err_comb) %% length(D) # Find min and max position values
    Xminmax <- matrix(approx(x = X, xout = Xpos_minmax, rule = 2)$y, ncol = 2) # Find associated D values (linear interpolation)
    x_err_comb <- 0.5 * ((Xminmax[,2] - Xminmax[,1]) %% X[length(X)]) # Find combined SD in D domain
    return(x_err_comb)
}