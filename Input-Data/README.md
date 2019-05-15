# Preparing Inputs

Here we store the scripts used to build the PHONEMeS inputs used for the analysis from data and their outputs. These inputs are generated simply by running the `R` scripts `prepareGMM_FC_UACC257_redundant_sites.R` and `prepareGMM_FC_A2058_redundant_sites.R` for each cell-lines respecitvely.

## Steps to build the inputs

Below we describe the main steps followed to build the PHONEMeS inputs.

+ We read the Mass-Spect data `inst/Endothelin _ShortTerm_UACC257_new180822_SiteLevel.csv` and `inst/Endothelin _ShortTerm_A2058_new_180822_SiteLevel.csv` containing the site level informations (fold changes and significance as described in the paper).
+ We `log2` transform the fold changes.
+ We assign a threshold value `threshQ=0.1` to the differential abundance.
+ Based on the data tables, we generate the GMM objects containing a list of matrices for each site. Each matrix contains different values and labels assigned to each measurement at all time-points as below:
+ The `Indiv` columns contain the scores we assign to each of the measurements by computing the `log2` value of the ratio between the significance of a measurement at a specific time-point (as found in the data tables) over the `threshQ` value. Sites with a `qValue < threshQ` on a specific time-point will be assigned a negative score, otherwise they will be assigned a positive score.
+ The `clus` columns contain information about significance status of the measurement. If `qValue < threshQ`, then to the data-point it will be assigned a Perturbation `P` status, otherwise the measurement will take a Control `C` status.
+ The `FCvCaPval` contain the `q` significance values for each measurement.
+ The `status` columns show whether a measurement has been observed to be higly regulated (`OK`) or not (`FP`). A site with a fold change higher than 1.5 is considered to be as highly regulated.
+ The `FCvC` column contains the log2 fold changes of EDN activated sites compared to control.

## Generating and storing the inputs

We build two separate data objects `dataGMM_UACC257.RData` and `dataGMM_A2058.RData` and we keep them stored on the local directory.
