# ShellChron 0.3.5

## 0.2.8
Added escape statements (on.exit()) to functions that change the WD
Updated description to refer to published manuscript with more details
on the methodology.

## 0.2.7
Removed characters that caused problems with building PDF manual
Minor spelling checks
Updated some formulae to make checking run more smoothly

## 0.3.5
Following peer-review, I have implemented options in ShellChron that enable the user to set the parameters of the SCEUA algorithm and disable the sinusoidal regression (sinreg).
Previous versions of ShellChron underestimated the seasonal δ18Oc range and returned sequences of constant δ18Oc values in individual modelling windows. These issues were caused by a combination of issues in the run_model and cumulative_day functions:
1. The boundaries on the growth rate sinusoid were set too liberally, causing growth rates to vary by several orders of magnitude within years and allowing growth rate variability over time to be nearly flat or vertical. The run_model function was modified to make the growth rate boundaries more restrictive, and dependent on the estimated initial value for the mean growth rate. The new boundaries still allow growth rate to vary by a factor 10 within a year, which should cover the natural variability of growth rates in δ18Oc archives.
2. Constant δ18Oc solutions in modelling windows were often a result of the number of years of δ18Oc that were modelled within a window. If the window contained more time than expected the samples outside the δ18Oc curve were assigned a constant value, resulting in a flat δ18Oc
profile. This issue also explained why the mean modeled δ18Oc values underestimated the true δ18Oc variability in the record.
3. Slight modifications were made to the cumulative_day function that aligns results from individual modelling windows to a common time axis. The previous version of ShellChron added too many or too little days to the age results in modelling windows in which either none or multiple year transitions occur. This resulted in sudden jumps in time in the age-distance relationship.