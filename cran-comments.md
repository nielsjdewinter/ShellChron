## Resubmission
This is a resubmission. In this version I have made the following changes (details in NEWS):
* Updated run_model function to prevent "flat lines" in d18Oc models
* Updated run_model and wrap_function function to allow more flexibility on model parameters
* Updated cumulative_day function to prevent glitches when adding time to windows to put them on cumulative timescale, causing jumps in time
* Removed appearances of the hyphen (U+2010) failing manual building from Rd files

## Test environments
* local Win 10 Pro install, R 4.0.2
* Ubuntu 16.04 (on travis ci) R 4.0.2
* OS X 11.3 (on travis ci) R 4.1.0

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## R-hub Builder results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
There are currently no downstream dependencies for this package