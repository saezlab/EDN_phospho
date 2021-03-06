# Polishing Networks

Here we store the scripts used for combining our network solutions as final results and label nodes for a nice visualization in Cytoscape.

## Steps to prepare the final results

+ We first combine all the solutions from the 100 iterations (stored in `../Analysis/resList_UACC257.RData` or `../Analysis/resList_A2058.RData`) into a single table.
+ We assign weights to each interaction in the combined solution based on how many times it has appeared across all the 100 models. Higher weight, means that the interaction has appeared more frequently and we are more confident about it's presence.
+ For each interaction we also assign a label based on which time-point it appears. Earlier time-points prevail over later ones when labeling the edges.
+ Here we also assign node attributes to each of the nodes based on what they represent: `D` for the EDNRB and `P` to measurements.
+ We then map the nodes present in the final solution from UniProt to the more common Gene ID's as based on the mapping table `MappingUniprot_Gene_names 16.46.06.csv`.
+ Finally we retain in the final solution the most confident interactions by removing the ones with assigned weight smaller than 20.
+ **The script `visualNetwork_polish.R` has to be executed for each cell-line by changing the commented lines in the script corresponding to the input data, network and output files.**

## Generating and storing the final results

Edge weights and time-point attributes together with node attributes are generated simply by running `assign_timepoint_bootstrap_AS_UACC257.R` or `assign_timepoint_bootstrap_AS_A2058.R` scripts. Mapping to gene ID's is done right after simply by executing the `visualNetwork_polish.R` script. Final combined and reduced network results are saved in the `../Results/` directory.
