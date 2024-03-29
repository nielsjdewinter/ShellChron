% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_model.r
\name{run_model}
\alias{run_model}
\title{Function that optimizes sinusoid parameters to fit d18O data}
\usage{
run_model(
  dat,
  dynwindow,
  transfer_function = "KimONeil97",
  d18Ow = 0,
  T_per = 365,
  G_per = 365,
  t_int = 1,
  t_maxtemp = 182.5,
  SCEUApar = c(1, 25, 10000, 5, 0.01, 0.01),
  sinfit = TRUE,
  MC = 1000,
  plot = TRUE
)
}
\arguments{
\item{dat}{Matrix containing the input data}

\item{dynwindow}{Information on the position and length of modeling
windows}

\item{transfer_function}{Transfer function used to convert d18Oc to temperature
data.}

\item{d18Ow}{Either a single value (constant d18Ow) or a vector of length
equal to the period in SST data (365 days by default) containing information
about seasonality in d18Ow. Defaults to constant d18Ow of 0 permille VSMOW
(the modern mean ocean value)}

\item{T_per}{Period of SST sinusoid (in days; default = 365)}

\item{G_per}{Period of growth rate sinusoid (in days; default = 365)}

\item{t_int}{Time interval (in days; default = 1)}

\item{t_maxtemp}{Timing of the warmest day of the year (in julian day;
default = 182.5, or May 26th halfway through the year)}

\item{SCEUApar}{Parameters for SCEUA optimization (iniflg, ngs, maxn, kstop
pcento, peps). For details, refer to Duan et al. (1992) in references}

\item{sinfit}{Apply sinusoidal fitting to guess initial parameters for SCEUA
optimization? \code{TRUE/FALSE}}

\item{MC}{Number of Monte Carlo simulations to apply for error propagation
Default = 1000}

\item{plot}{Should results of modeling be plotted? \code{TRUE/FALSE}}
}
\value{
A list containing the \code{resultarray} which contains the full
result of all simulations on each data window and the \code{parmat} listing
all optimized growth rate and SST parameters used to model d18O in each data
window
}
\description{
The second core function of the ShellChron growth model. Loops
through all data windows and uses the \code{growth_model} function
to create d18O series that match the input data. This step is
iterated and optimized (minimizing the Sum of Squared Residuals)
through the SCEUA algorithm (by Duan et al., 1992) which finds
the optimal input parameters to the growth rate and Sea Surface
Temperature (SST) sinusoids to simulate d18O data.
}
\examples{
# Create dummy input data column by column
dat <- as.data.frame(seq(1000, 40000, 1000))
colnames(dat) <- "D"
dat$d18Oc <- sin((2 * pi * (seq(1, 40, 1) - 8 + 7 / 4)) / 7)
dat$YEARMARKER <- c(0, rep(c(0, 0, 0, 0, 0, 0, 1), 5), 0, 0, 0, 0)
dat$D_err <- rep(100, 40)
dat$d18Oc_err <- rep(0.1, 40)
# Create dummy dynwindow data
dynwindow <- as.data.frame(seq(1, 29, 2))
colnames(dynwindow) <- "x"
dynwindow$y <- rep(12, 15)
# Run model function
\donttest{resultlist <- run_model(dat = dat,
    dynwindow = dynwindow,
    transfer_function = "KimONeil97",
    d18Ow = 0,
    T_per = 365,
    G_per = 365,
    t_int = 1,
    t_maxtemp = 182.5,
    SCEUApar = c(1, 25, 10000, 5, 0.01, 0.01),
    sinfit = TRUE,
    MC = 1000,
    plot = FALSE)}
}
\references{
package dependencies: ggplot2 3.2.1; rtop 0.5.14
Function dependencies: sinreg, d18O_model, growth_model

\doi{10.1029/91WR02985}
}
\seealso{
Duan, Qingyun, Soroosh Sorooshian, and Vijai Gupta. "Effective and
efficient global optimization for conceptual rainfall runoff models." Water
resources research 28.4 (1992): 1015-1031. https://doi.org/10.1029/91WR02985
}
