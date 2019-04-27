# Building Background network

PHONEMeS is a method that trains a prior knowledge network to the data. To adapt the method to GCPR signalling the background network (representing our prior knowledge), was build as follows: First, the background kinase-to-substrate (K-S) network derived from [Omnipath](http://omnipathdb.org/) was complemented with a subset of all directed and signed PPIs from [Omnipath](http://omnipathdb.org/) associated with GPCR signalling in Reactome (Pathway R-HSA-372790). Closely related heterotrimeric G-Protein subunits and kinase isoforms were grouped in the prior knowledge database.

## Steps to build the background network

Below we describe the main resources used to build our background network.

+ We start by adding all the Omnipath K-S interactions as stored in `inst/ptms.txt` file.
+ We identify the protein associated with GPCR signalling in Reactome ([Pathway R-HSA-372790]((http://software.broadinstitute.org/gsea/msigdb/cards/REACTOME_GPCR_DOWNSTREAM_SIGNALING))). This list of proteins is the one stored in `inst/geneset.txt` file.
+ From the list of protein interactions of Omnipath (`inst/Omnipath_interactions.txt`) we select only those that are signed and directed and which are involving the proteins associated to GPCR signalling.
+ We can add or remove certain Omnipath interactions based on the reasons explaied in the `inst/IntMod.csv` and `inst/PhosMod.csv` files.
+ We group some isoform proteins based on how similar are they structurally and functionally. Grouped isoforms are described in `inst/Group.csv` file.

## Generating and storing the background network

Our background network is generated simply by executing the `prepareBN.R` script and which applies all the steps described above. The output of the analysis is the `interactionsALL.RData` object which contains the list of all 15097 interactions used for the training of the models for both the cell lines.

## Side note

We have used the Omnipath version of *March 2018* to build our background network.
