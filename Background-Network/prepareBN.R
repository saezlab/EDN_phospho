#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga & Alexander Schaefer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmai.com
#
##############################################################################
# 12:58 24/03/2018
# This script is used to produce the background network which is then used to
# train the PHONEMeS networks.

library(readr)

# Loading Omnipath interactions and ptms lists as of March 2018
Omnipath_interactions <- read_delim("inst/Omnipath_interactions.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # Loading omnipath interactions
ptms <- read_delim("inst/ptms.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # Loading the ptms tables

allSpecies <- unique(c(Omnipath_interactions$source, Omnipath_interactions$target, ptms$enzyme, ptms$substrate))
write.table(x = allSpecies, file = "inst/toUniprot.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE) # file to map the uniprot ID's: https://www.uniprot.org/

uniprot <- read_delim("inst/uniprot.tab", "\t", escape_double = FALSE, trim_ws = TRUE) # Uniprot table after the mapping

######
# Adding PTMS interactions
bn <- matrix(data = , nrow = nrow(ptms), ncol = 8) # Initializing empty BN
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

idx <- which(ptms$modification%in%c("phosphorylation")) # we only keep the phosphorylations
ptms <- ptms[idx, ]

cnt <- 1
for(i in 1:nrow(ptms)){
  
  bn[i, 1] <- as.character(ptms$substrate[i])
  bn[i, 2] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 1])[1]]
  bn[i, 3] <- as.character(ptms$enzyme[i])
  bn[i, 4] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 3])[1]]
  bn[i, 5] <- ptms$residue_type[i]
  bn[i, 6] <- as.character(ptms$residue_offset[i])
  bn[i, 7] <- paste0("e", cnt)
  cnt <- cnt + 1
  bn[i, 8] <- paste0(bn[i, 2], "_", bn[i, 5], bn[i, 6])
  
}

BN1 <- bn
nrow(BN1) # adding 15623 ptms interactions

######
# Adding Omnipath interactions
sum1 <- Omnipath_interactions$is_stimulation
sum2 <- Omnipath_interactions$is_inhibition
sum <- sum1 + sum2
idx <- which(sum>=1)
omnipath <- Omnipath_interactions[idx, ] #removing PPI's and keeping only the directed and signed interactions

geneset <- read_csv("inst/geneset.txt", col_names = FALSE) # Loading the Genes of the MSigDB database for the GPCR downstream signalling: http://software.broadinstitute.org/gsea/msigdb/cards/REACTOME_GPCR_DOWNSTREAM_SIGNALING

genes2uniprot <- read_delim("inst/genes2uniprot.tab", "\t", escape_double = FALSE, trim_ws = TRUE) # Loading a mapping table of genes to uniprot ID's from uniprot

uGenes <- unique(genes2uniprot$`Gene names`)

uGenes <- unique(unlist(lapply(strsplit(x = uGenes, split = " "), `[[`, 1)))
mapGenes <- matrix(data = , nrow = length(uGenes), ncol = 3)
mapGenes[, 1] <- uGenes
for(i in 1:length(uGenes)){
  
  geneID <- paste0(uGenes[i], "_HUMAN")
  if(geneID%in%genes2uniprot$`Entry name`){
    mapGenes[i, 2] <- geneID
    mapGenes[i, 3] <- uniprot$Entry[which(uniprot$`Entry name`==geneID)[1]]
  } else {
    if(length(which(genes2uniprot$`Gene names`==uGenes[i])) > 0){
      mapGenes[i, 2] <- genes2uniprot$`Entry name`[which(genes2uniprot$`Gene names`==uGenes[i])[length(which(genes2uniprot$`Gene names`==uGenes[i]))]]
      mapGenes[i, 3] <- uniprot$Entry[which(uniprot$`Entry name`==mapGenes[i, 2])[1]] 
    }
  }
  
}

mapGenes <- mapGenes[complete.cases(mapGenes), ]

# now selecting from Omnipath pathways only those interactions involving both of the mapped genes of the GPCR downstream signalling
idx1 <- which(omnipath$source%in%mapGenes[, 3])
idx2 <- which(omnipath$target%in%mapGenes[, 3])
idx <- intersect(x = idx1, y = idx2)

bn <- matrix(data = , nrow = length(idx), ncol = 8) # Initializing empty BN
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

for(i in 1:nrow(bn)){
  
  bn[i, 1] <- as.character(omnipath$target[idx[i]])
  bn[i, 2] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 1])[1]]
  bn[i, 3] <- as.character(omnipath$source[idx[i]])
  bn[i, 4] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 3])[1]]
  bn[i, 5] <- "R"
  bn[i, 6] <- "1"
  bn[i, 7] <- paste0("e", cnt)
  cnt <- cnt + 1
  bn[i, 8] <- paste0(bn[i, 2], "_R1")
  
}

BN2 <- bn
nrow(bn) # adding 270 omnipath interactions


#####
# adding RAS interactions from Omnipath

ras <- c(uniprot$Entry[which(uniprot$`Entry name`=="RASK_HUMAN")], 
         uniprot$Entry[which(uniprot$`Entry name`=="RASH_HUMAN")],
         uniprot$Entry[which(uniprot$`Entry name`=="RASN_HUMAN")])

idx <- unique(c(which(omnipath$source%in%ras), which(omnipath$target%in%ras)))

bn <- matrix(data = , nrow = length(idx), ncol = 8) # Initializing empty BN
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

for(i in 1:length(idx)){
  
  bn[i, 1] <- omnipath$target[idx[i]]
  bn[i, 2] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 1])[1]]
  bn[i, 3] <- omnipath$source[idx[i]]
  bn[i, 4] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 3])[1]]
  bn[i, 5] <- "R"
  bn[i, 6] <- "1"
  bn[i, 7] <- paste0("e", cnt)
  cnt <- cnt + 1
  bn[i, 8] <- paste0(bn[i, 2], "_R1")
  
}

nrow(bn)

BN2 <- rbind(BN2, bn)

######
# Suggested Modifications - 1 (ADD & REMOVE Interactions)
IntMod <- read_csv("inst/IntMod.csv")
idx <- which(IntMod$Modification=="REMOVE")
temp <- bn
ii <- c()
for(i in 1:length(idx)){
  
  ii <- c(ii, intersect(which(BN2[, 3]==IntMod$source[idx[i]]), which(BN2[, 1]==IntMod$target[idx[i]])))
  
}

BN2 <- BN2[-ii, ]
length(ii) # 5 interaction removed

idx <- which(IntMod$Modification=="ADD")

bn <- matrix(data = , nrow = length(idx), ncol = 8) # Initializing empty BN
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")
for(i in 1:length(idx)){
  
  bn[i, 1] <- as.character(IntMod$target[idx[i]])
  bn[i, 2] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 1])[1]]
  bn[i, 3] <- as.character(IntMod$source[idx[i]])
  bn[i, 4] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 3])[1]]
  bn[i, 5] <- "R"
  bn[i, 6] <- "1"
  bn[i, 7] <- paste0("e", cnt)
  cnt <- cnt+1
  bn[i, 8] <- paste0(bn[i, 2], "_R1")
  
}

BN3 <- bn
nrow(bn) # 44 more interactions added

######
# Suggested Modifications - 2 (ADD & REMOVE Phos)
PhosMod <- read_csv("inst/PhosMod.csv")
PhosRm<-subset(PhosMod, Mod == "REMOVE")

ii <- c()
for(i in 1:nrow(PhosRm)){
  
  idx1 <- intersect(x = which(BN1[, 3]==PhosRm$Kinase[i]), y = which(BN1[, 1]==PhosRm$Substrate[i]))
  idx2 <- intersect(x = which(BN1[, 5]==PhosRm$Site[i]), y = which(BN1[, 6]==PhosRm$Residue[i]))
  
  ii <- c(ii, intersect(x = idx1, y = idx2))
  
}

BN1 <- BN1[-ii, ]
length(ii) # removing 73 k-s interactions

PhosAdd<-subset(PhosMod, Mod == "ADD")

bn <- matrix(data = , nrow = nrow(PhosAdd), ncol = 8) # Initializing empty BN
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

cnt <- 1
for(i in 1:nrow(PhosAdd)){
  
  bn[i, 1] <- as.character(PhosAdd$Substrate[i])
  bn[i, 2] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 1])[1]]
  bn[i, 3] <- as.character(PhosAdd$Kinase[i])
  bn[i, 4] <- uniprot$`Entry name`[which(uniprot$Entry==bn[i, 3])[1]]
  bn[i, 5] <- PhosAdd$Site[i]
  bn[i, 6] <- as.character(PhosAdd$Residue[i])
  bn[i, 7] <- paste0("e", cnt)
  cnt <- cnt + 1
  bn[i, 8] <- paste0(bn[i, 2], "_", bn[i, 5], bn[i, 6])
  
}

BN4<-bn

######
# Suggested Modifications - 3 (Adding EDNRB interactions)
Group <- read_csv("inst/Group.csv")
Gprot<-Group[1:34,]
g <- unique(Gprot$`Entry name`)

bn <- matrix(data = , nrow = length(g), ncol = 8) # adding EDNRB~Gprotein interactions
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

i<-1
cnt<-1
for(i in 1:length(g)){
  
  bn[i, 1] <- Gprot$Entry[i]
  bn[i, 2] <- g[i]
  bn[i, 3] <- "P24530"
  bn[i, 4] <- "EDNRB_HUMAN"
  bn[i, 5] <- "R"
  bn[i, 6] <- "1"
  bn[i, 7] <- paste0("e", cnt)
  cnt <- cnt + 1
  bn[i, 8] <- paste0(g[i], "_R1")
  
}

BN5 <- bn
nrow(bn) # 34 interactions added


BN <- rbind(BN1, BN2, BN3, BN4, BN5)


nrow(BN)
BN <- BN[complete.cases(BN), ]
nrow(BN)


## Grouping the RAS isoforms

BN[, 2] <- gsub(pattern = "RASH_HUMAN", replacement = "RAS_HUMAN", x = BN[, 2])
BN[, 4] <- gsub(pattern = "RASH_HUMAN", replacement = "RAS_HUMAN", x = BN[, 4])
BN[, 8] <- gsub(pattern = "RASH_HUMAN", replacement = "RAS_HUMAN", x = BN[, 8])

BN[, 2] <- gsub(pattern = "RASK_HUMAN", replacement = "RAS_HUMAN", x = BN[, 2])
BN[, 4] <- gsub(pattern = "RASK_HUMAN", replacement = "RAS_HUMAN", x = BN[, 4])
BN[, 8] <- gsub(pattern = "RASK_HUMAN", replacement = "RAS_HUMAN", x = BN[, 8])

BN[, 2] <- gsub(pattern = "RASN_HUMAN", replacement = "RAS_HUMAN", x = BN[, 2])
BN[, 4] <- gsub(pattern = "RASN_HUMAN", replacement = "RAS_HUMAN", x = BN[, 4])
BN[, 8] <- gsub(pattern = "RASN_HUMAN", replacement = "RAS_HUMAN", x = BN[, 8])


# Now grouping the isoforms
for(i in 1:nrow(Group)){
  
  idx <- which(BN[, 2]==as.character(Group$`Entry name`[i]))
  if(length(idx) > 0){
    BN[idx, 2] <- Group$Group[i]
    BN[idx, 8] <- paste0(BN[idx, 2], "_", BN[idx,5], BN[idx,6])
  }
  idx <- which(BN[, 4]==as.character(Group$`Entry name`[i]))
  if(length(idx) > 0){
    BN[idx, 4] <- Group$Group[i]
  }
  
}


nrow(BN) # 15610 in the final bacground network

# Now building the Background Network as a data-frame
allD <- as.data.frame(BN)
allD$S.AC <- as.character(allD$S.AC)
allD$S.ID <- as.character(allD$S.ID)
allD$K.AC <- as.character(allD$K.AC)
allD$K.ID <- as.character(allD$K.ID)
allD$res <- as.character(allD$res)
allD$pos <- as.character(allD$pos)
allD$SID <- paste0("e", 1:nrow(allD))
allD$S.cc <- as.character(allD$S.cc)

# Saving the background-network as data-frame
save(allD, file = "interactionsALL.RData")

