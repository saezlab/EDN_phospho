# Analysis

This directory contains scripts necessary to perform the PHONEMeS analysis. PHONEMeS has been reformulated as an Integer Linear Programming (ILP) problem and this new version is the one used for this study. **For performing th ILP analysis, users must first obtain a [CPLEX license](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_%7C470%7C135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB) (free for academic purposes) and store on this working directory the cplex executable file.**

## Steps to run the analysis

+ First we load our generated background network and data inputs from the previous pre-processing steps.
+ We then load some of the functions used for the time-point PHONEMeS analysis and which are stored in the `/Public` directory.
+ We specify the names of our conditions as we have them on the `GMM` and `GMM.all` generated objects and in this case, for each condition where we have one time-point.
+ We run the `PHONEMeS_dt` analysis 100 times where for each iteration we downsample the data. In this case, the analysis is performed by generating an ILP problem which is then read and optimized by the `cplex` solver. All the functions used to generate the problem, optimize it and retrieve the results are called from the `/Public` directory.
+ The 100 time-specific models are generated, combined and saved as a list

## Generating and storing the analysis results

All it takes to perform the anlysis of the networks is to simply run the `executionScript.R`. The analysis might take up to some minutes until all the iterations are finished. The output of the analysis are saved locally as `resList_UACC257.RData` or `resList_A2058.RData`. These data objects contain the resulting network obtained from the combination of each time-point and iteration.
