#' Function that identifies peaks in a dataset
#' 
#' Developed by William A. Huber
#' @param x Vector of \code{x} values of input data
#' @param y Vector of \code{y} values of input data
#' @param w Window size for smoothing data
#' @return A vector listing the standard deviations of propagated errors 
#' propagated on all \code{X} values.
#' @seealso https://rpubs.com/mengxu/peak_detection
#' @references package dependencies: zoo 1.8.7
#' @examples
#' # Create dummy periodic data
#' x <- seq(1, 100, 1)
#' y <- sin((2 * pi * (seq(1, 100, 1) - 8 + 20 / 4)) / 20)
#' # Run peakID function
#' result <- PeakID(x, y, w = 20)
#' @export
PeakID <- function(x, y, w=1, ...) { # Algorhythm for peak identification by William A. Huber (https://rpubs.com/mengxu/peak_detection)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- zoo::rollapply(zoo::zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}