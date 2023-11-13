# MODGEE model

## About the code

This code implements the 'matrix of detections gives extinction estimates' or MODGEE model to estimate the number of plant, bird, and butterfly extinctions that have occurred in Singapore, including species that went extinct before they could be recorded, i.e., "dark extinctions".

## About the manuscript

Chisholm, R.A. (2023) Two centuries of biodiversity discovery and loss in Singapore, *PNAS*

## Abstract

(coming soon)

## Notes for running the code

- The "SEUX" folder has a script `SEUX_Singapore.R` that runs all the SEUX analyses. 
- The "MODGEE" folder has a script `MODGEE_Singapore.R` that runs the MODGEE analyses. One must modify a line of code near the beginning to indicate which taxonomic group to run.
- The `get_CI_estimate()` function has been modified from the orignal function in the "seux" library to return the complete data from all the repeats.
- For the purposes of quickly demonstrating how the estimates are obtained, the code prints only the estimate from the main model and therefore will not match the published results. In order to obtain the full published results, the mean was calculated over 10,000 repeats from the confidence interval estimator. Please contact us if you need assistance reproducing our main results.

## Related projects

More information on the SEUX method, which requires as input only first- and last-record dates, can be found here: https://github.com/nadiahpk/seux

