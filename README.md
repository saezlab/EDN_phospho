# EDN_phospho

Elucidating essential kinases of endothelin signalling by logic modelling of phosphoproteomics data

## Purpose

Here I provide documentation and the analysis scripts used for generating the time-resolved network models of EDN signalling in UACC257 and
A2058 cell-lines.

## Organisation of the repository

+ [Preparing inputs](https://github.com/saezlab/EDN_phospho/tree/master/Input-Data): Here we prepare the PHONEMeS inputs objects from the processed data.
+ [Building background networks](https://github.com/saezlab/EDN_phospho/tree/master/Background-Network): Here we prepare the background network which is used for the training EDN signalling to data.
+ [Analysis](): Here we perform the PHONEMeS time-point analysis to train our models.
+ [Results](): Here we store the analysis results.

## Requirements

All network modelling steps were performed in R v3.4 and visualised using Cytoscape v3.3.

For running the PHONEMeS analysis, the user must obtain the cplex [license](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_%7C470%7C135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB)
from IBM, which is for free for academic purposes. Once the user obtains the license, they should store the cplex executable file to the 
working directory where PHONEMeS is running (here: on the [Analysis]() directory).

Other prerequisites include downloading and installing the following `R` package dependencies:

+ [PHONEMeS](https://saezlab.github.io/PHONEMeS/)
+ [igraph](https://igraph.org/r/)
+ [BioNet](https://bioconductor.org/packages/release/bioc/html/BioNet.html)
+ [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)

The Analysis was run under MAC OS X El Captain machine. Please note that the results might slightly change when running the anlysis from
machines with different OS.
