#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmai.com
#
##############################################################################
# 12:59 23/03/2018
# Asigning nodes attributes for better visualization of PHONEMeS networks

load(file = "../Input-Data/dataGMM_UACC257.RData") # for uacc257
#load(file = "../Input-Data/dataGMM_A2058.RData") # for a2058

nodesMatrix <- matrix(data = , nrow = (length(GMM.ID$S.cc)+1), ncol = 2)
nodesMatrix[1, 1] <- "EDNRB_HUMAN"
nodesMatrix[1, 2] <- "D"
nodesMatrix[2:(nrow(nodesMatrix)), 1] <- as.character(GMM.ID$S.cc)
nodesMatrix[2:(nrow(nodesMatrix)), 2] <- "P"

colnames(nodesMatrix) <- c("Species", "nodesP")
write.table(x = nodesMatrix, file = "inst/nodesAttributes.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

network<-"inst/tpBootstrapAnalysis_UACC257.txt" # for uacc257
#network<-"inst/tpBootstrapAnalysis_A2058.txt" # for a2058
att<-"inst/nodesAttributes.txt"

SIF<-as.matrix(read.table(network, header =TRUE))
SIF<-gsub("_HUMAN","",SIF)
x<-SIF[grep("_",SIF[,3]),3]
x<-x[!grepl("GNA", x)] 
x<-cbind(x,rep("S",length(x)))
nodes<-as.matrix(read.table(att, header =TRUE))
nodes<-gsub("_HUMAN","",nodes)
nodes2<-rbind(nodes,x)

nodes3<-cbind(unique(nodes2[,1]),"")
nodes3[,2]<-nodes2[,2][match(nodes3[,1], nodes2[,1])]
colnames(nodes3)<-colnames(nodes2)

##
# Mapping network nodes
library(readr)

mapping = read_csv("inst/MappingUniprot_Gene_names 16.46.06.csv")

sif = SIF

for(ii in 1:nrow(x = SIF)){
  
  idx = which(mapping$`name after script (Uniprot)`==SIF[ii, 1])
  if(length(idx)>0){
    sif[ii, 1] = mapping$`name in paper (Gene names)`[idx[length(idx)]]
  }
  
  idx = which(mapping$`name after script (Uniprot)`==SIF[ii, 3])
  if(length(idx)>0){
    sif[ii, 3] = mapping$`name in paper (Gene names)`[idx[length(idx)]]
  }
  
}

nodes_attrib = nodes3
for(ii in 1:nrow(mapping)){
  
  idx = which(nodes_attrib[, 1]==mapping$`name after script (Uniprot)`[ii])
  if(length(idx) > 0){
    nodes_attrib[idx, 1] = mapping$`name in paper (Gene names)`[ii]
  }
}

write.table(sif, file = "../Results/sif_uacc257.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(sif[which(as.numeric(sif[, 4])>20), ], file = "../Results/sif_uacc257_red.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(nodes_attrib, file = "../Results/attrib_uacc257.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

#write.table(sif, file = "../Results/sif_a2058.txt", quote = FALSE, row.names = FALSE, sep = "\t")
#write.table(sif[which(as.numeric(sif[, 4])>20), ], file = "../Results/sif_a2058_red.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(nodes_attrib, file = "../Results/attrib_a2058.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

