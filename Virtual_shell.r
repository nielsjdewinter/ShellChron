#' Virtual d18O data of a 6 year shell record
#'
#' A dataset containing d18O values for 6 seasons of d18O
#' data with stochastic noise to simulate random weather
#' variability on the seasonality.
#'
#' @format A matrix with 80 rows and 5 variables:
#' \describe{
#'   \item{D}{depth, in micrometers}
#'   \item{d18Oc}{stable oxygen isotope ratio of shell carbonate,
#'    in permille VPDB}
#'   \item{YEARMARKER}{marks year transitions with a "1"}
#'   \item{D_err}{depth uncertainty, in micrometers}
#'   \item{d18Oc_err}{d18O uncertainty, in permille VPDB}
#' }
#' @source \url{https://cp.copernicus.org/preprints/cp-2020-118/}
"Virtual_shell"