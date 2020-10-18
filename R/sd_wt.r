#' Function to calculate weighted standard deviation
#' 
#' Calculates the standard deviation of a weighted sample set
#' while propagating sample weights through the calculation.
#' @param x Vector containing the values in the set
#' @param w Vector containing the weights to each value (in 
#' the same order as \code{x}
#' the optimized growth model
#' @param na.rm Should NA values be removed from the set prior
#' to calculation? \code{TRUE/FALSE}
#' @return The standard deviation of the weighted set of \code{x} values
#' @examples
#' # Create dummy data
#' x <- seq(1, 10, 0.5)
#' w <- c(seq(0.1, 1, 0.1), seq(0.9, 0.1, -0.1))
#' SDw <- sd_wt(x, w, na.rm = TRUE) # Run the function
#' @export
sd_wt<-function(x, w, na.rm = FALSE){ # Formula for weighted standard deviation
    if(na.rm == TRUE){ # Remove NA containing x/w pairs
        x<-x[!(is.na(x) | is.na(w))]
        w<-w[!(is.na(x) | is.na(w))]
    }
    mean.wt <- mean(x * w, na.rm = TRUE) / mean(w) # Calculate weighted mean
    stdev <- sqrt(sum(w * (x - mean.wt) ^ 2) / ((length(w) - 1) / length(w) * sum(w)))
    return(stdev)
}