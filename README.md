# MODGEE and SEUX estimates of Singapore extinctions

## About the code

This code implements two statical methods for estimating the number of species that have occurred in Singapore, including species that went extinct before they could be recorded, i.e., "dark extinctions": the 'matrix of detections gives extinction estimates' or MODGEE model, and the SEUX model.

## About the manuscript

Chisholm, R.A. _et al._ (2023) Two centuries of biodiversity discovery and loss in Singapore, *PNAS*

## Abstract

There is an urgent need for reliable data on the impacts of deforestation on tropical biodiversity. The city-state of Singapore has one of the most detailed biodiversity records in the tropics, dating back to the turn of the 19th century. In 1819, Singapore was almost entirely covered in primary forest, but this has since been largely cleared. We compiled more than 200 y of records for 10 major taxonomic groups in Singapore (>50,000 individual records; >3,000 species), and we estimated extinction rates using recently developed and novel statistical models that account for “dark extinctions,” i.e., extinctions of undiscovered species. The estimated overall extinction rate was 37% (95% CI [31 to 42%]). Extrapolating our Singapore observations to a future business-as-usual deforestation scenario for Southeast Asia suggests that 18% (95% CI [16 to 22%]) of species will be lost regionally by 2100. Our extinction estimates for Singapore and Southeast Asia are a factor of two lower than previous estimates that also attempted to account for dark extinctions. However, we caution that particular groups such as large mammals, forest-dependent birds, orchids, and butterflies are disproportionately vulnerable.

## Notes for running the code

- The "SEUX" folder has a script `SEUX_Singapore.R` that runs all the SEUX analyses. 
- The "MODGEE" folder has a script `MODGEE_Singapore.R` that runs the MODGEE analyses. One must modify a line of code near the beginning to indicate which taxonomic group to run.
- The `get_CI_estimate()` function has been modified from the orignal function in the "seux" library to return the complete data from all the repeats.
- For the purposes of quickly demonstrating how the estimates are obtained, the code prints just the estimate of fraction of extinct species from the main model, which will not match the published results exactly. To generate the published results exactly, users should instead take the mean over 10,000 repeats from the confidence interval estimator.

## Related projects

More information on the SEUX method, which requires as input only first- and last-record dates, can be found here: https://github.com/nadiahpk/seux

