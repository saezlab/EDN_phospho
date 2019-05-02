#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga & Alexander Schaefer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmail.com
#
##############################################################################
# 14:01 05/02/2019
# This script is used to produce the GMM objects which are the used as an input
# for PHONEMeS functions and analysis -- Endothelin _ShortTerm_A2058_new_180822_SiteLevel.csv

library(readr)
Endothelin_ShortTerm_A2058 <- read_csv("inst/Endothelin _ShortTerm_A2058_new_180822_SiteLevel.csv")

wtFC <- Endothelin_ShortTerm_A2058[, c("2minFC", "10minFC", "20minFC", "30minFC", "60minFC", "90minFC")]

wtQ <- Endothelin_ShortTerm_A2058[, c("2min q", "10min q", "20min q", "30min q", "60min q", "90min q")]

rNames <- paste0(Endothelin_ShortTerm_A2058$"Uniprot accession", "_", Endothelin_ShortTerm_A2058$"Site", "_", Endothelin_ShortTerm_A2058$"n PO4")

rownames(wtFC) <- rNames
rownames(wtQ) <- rNames

wtFC[which(wtFC=="#DIV/0!", arr.ind = TRUE)] <- NA

###
# Preparing mapping table

data.IDmap <- matrix(, nrow = nrow(wtFC), ncol = 3)
colnames(data.IDmap) <- c("dataID", "UPID", "S.cc")
data.IDmap[, 1] <- rNames
data.IDmap[, 2] <- Endothelin_ShortTerm_A2058$"Uniprot accession"
data.IDmap[, 3] <- paste0(Endothelin_ShortTerm_A2058$"Uniprot accession", "_", Endothelin_ShortTerm_A2058$Site)

###
# Preparing values for GMM object

numData.Avg <- matrix(, nrow = nrow(wtFC), ncol = 6)
numData.Avg[, 1] <- log(x = as.numeric(as.matrix(wtFC)[, 1]), base = 2)
numData.Avg[, 2] <- log(x = as.numeric(as.matrix(wtFC)[, 2]), base = 2)
numData.Avg[, 3] <- log(x = as.numeric(as.matrix(wtFC)[, 3]), base = 2)
numData.Avg[, 4] <- log(x = as.numeric(as.matrix(wtFC)[, 4]), base = 2)
numData.Avg[, 5] <- log(x = as.numeric(as.matrix(wtFC)[, 5]), base = 2)
numData.Avg[, 6] <- log(x = as.numeric(as.matrix(wtFC)[, 6]), base = 2)
rownames(numData.Avg) <- rownames(wtFC)

allNumData <- c(as.numeric(numData.Avg[, 1]), as.numeric(numData.Avg[, 2]), as.numeric(numData.Avg[, 3]), as.numeric(numData.Avg[, 4]), as.numeric(numData.Avg[, 5]), as.numeric(numData.Avg[, 6]))
threshRHS <- median(allNumData, na.rm = TRUE)+2.5*mad(allNumData, na.rm = TRUE)
threshLHS <- median(allNumData, na.rm = TRUE)-2.5*mad(allNumData, na.rm = TRUE)
min<-apply(numData.Avg,MARGIN = 1, FUN = min, na.rm = TRUE)
max<-apply(numData.Avg,MARGIN = 1, FUN = max, na.rm = TRUE)

threshQ <- 0.1
##
tp1.p <- as.numeric(as.matrix(wtQ)[, 1])          # assigning corresponding q-value to each measurement
tp1.p[which(is.na(tp1.p))] <- 1                   # setting missing q-values to 1
tp1.lo <- rep(0, nrow(wtFC))                      # initializing scoring vector to be assigned to each measurement
names(tp1.lo) <- rownames(wtFC)
tp1.c <- rep("C", length(tp1.p))                  # initializing status vector to be assigned to each measurement (Control - C/Perturbed P)
tp1.c[which(tp1.p < threshQ)] <- "P"              # assigning the Perturbed - P status to measurements with q-value < threshQ
tp1.lo <- log(x = tp1.p/threshQ, base = 2)        # assigning scores to each measurement (negative score for P nodes and positive scores for C nodes)
table(tp1.c)
tp1.s <- rep("OK", length(tp1.p))                 # initializing regulation status
names(tp1.s) <- rownames(numData.Avg)             # giving FP status to measurements with lower fc than 1.5 and OK status to measurement with higher fc than 1.5
tp1.s[which(min > -0.58)] <- "FP"
tp1.s[which(max > 0.58)] <- "OK"
# tp1.s[which(abs(numData.Avg[,1]) < 1)] <- "FP"
for(i in 1:length(tp1.p)){
  if(tp1.c[i]=="C"){
    tp1.s[i] <- "OK"
  }
}
table(tp1.s)

##
tp2.p <- as.numeric(as.matrix(wtQ)[, 2])
tp2.p[which(is.na(tp2.p))] <- 1
tp2.lo <- rep(0, nrow(wtFC))
names(tp2.lo) <- rownames(wtFC)
tp2.c <- rep("C", length(tp2.p))
tp2.c[which(tp2.p < threshQ)] <- "P"
tp2.lo <- log(x = tp2.p/threshQ, base = 2)
table(tp2.c)
tp2.s <- rep("OK", length(tp2.p))
names(tp2.s) <- rownames(numData.Avg)
tp2.s[which(min > -0.58)] <- "FP"
tp2.s[which(max > 0.58)] <- "OK"
# tp2.s[which(abs(numData.Avg[,2]) < 1)] <- "FP"
for(i in 1:length(tp2.p)){
  if(tp2.c[i]=="C"){
    tp2.s[i] <- "OK"
  }
}
table(tp2.s)

##
tp3.p <- as.numeric(as.matrix(wtQ)[, 3])
tp3.p[which(is.na(tp3.p))] <- 1
tp3.lo <- rep(0, nrow(wtFC))
names(tp3.lo) <- rownames(wtFC)
tp3.c <- rep("C", length(tp3.p))
tp3.c[which(tp3.p < threshQ)] <- "P"
tp3.lo <- log(x = tp3.p/threshQ, base = 2)
table(tp3.c)
tp3.s <- rep("OK", length(tp3.p))
names(tp3.s) <- rownames(numData.Avg)
tp3.s[which(min > -0.58)] <- "FP"
tp3.s[which(max > 0.58)] <- "OK"
# tp3.s[which(abs(numData.Avg[,3]) < 1)] <- "FP"
for(i in 1:length(tp3.p)){
  if(tp3.c[i]=="C"){
    tp3.s[i] <- "OK"
  }
}
table(tp3.s)

##
tp4.p <- as.numeric(as.matrix(wtQ)[, 4])
tp4.p[which(is.na(tp4.p))] <- 1
tp4.lo <- rep(0, nrow(wtFC))
names(tp4.lo) <- rownames(wtFC)
tp4.c <- rep("C", length(tp4.p))
tp4.c[which(tp4.p < threshQ)] <- "P"
tp4.lo <- log(x = tp4.p/threshQ, base = 2)
table(tp4.c)
tp4.s <- rep("OK", length(tp4.p))
names(tp4.s) <- rownames(numData.Avg)
tp4.s[which(min > -0.58)] <- "FP"
tp4.s[which(max > 0.58)] <- "OK"
# tp4.s[which(abs(numData.Avg[,4]) < 1)] <- "FP"
for(i in 1:length(tp4.p)){
  if(tp4.c[i]=="C"){
    tp4.s[i] <- "OK"
  }
}
table(tp4.s)

##
tp5.p <- as.numeric(as.matrix(wtQ)[, 5])
tp5.p[which(is.na(tp5.p))] <- 1
tp5.lo <- rep(0, nrow(wtFC))
names(tp5.lo) <- rownames(wtFC)
tp5.c <- rep("C", length(tp5.p))
tp5.c[which(tp5.p < threshQ)] <- "P"
tp5.lo <- log(x = tp5.p/threshQ, base = 2)
table(tp5.c)
tp5.s <- rep("OK", length(tp5.p))
names(tp5.s) <- rownames(numData.Avg)
tp5.s[which(min > -0.58)] <- "FP"
tp5.s[which(max > 0.58)] <- "OK"
# tp5.s[which(abs(numData.Avg[,5]) < 1)] <- "FP"
for(i in 1:length(tp5.p)){
  if(tp5.c[i]=="C"){
    tp5.s[i] <- "OK"
  }
}
table(tp5.s)

##
tp6.p <- as.numeric(as.matrix(wtQ)[, 6])
tp6.p[which(is.na(tp6.p))] <- 1
tp6.lo <- rep(0, nrow(wtFC))
names(tp6.lo) <- rownames(wtFC)
tp6.c <- rep("C", length(tp6.p))
tp6.c[which(tp6.p < threshQ)] <- "P"
tp6.lo <- log(x = tp6.p/threshQ, base = 2)
table(tp6.c)
tp6.s <- rep("OK", length(tp6.p))
names(tp6.s) <- rownames(numData.Avg)
tp6.s[which(min > -0.58)] <- "FP"
tp6.s[which(max > 0.58)] <- "OK"
# tp6.s[which(abs(numData.Avg[,6]) < 1)] <- "FP"
for(i in 1:length(tp6.p)){
  if(tp6.c[i]=="C"){
    tp6.s[i] <- "OK"
  }
}
table(tp6.s)

###
# Preparing GMM objects
GMM<-vector("list", dim(wtFC)[1])
names(GMM)<-rownames(wtFC)
for(i in 1:length(GMM)){
  GMM[[i]]<-rbind(
    c(as.character(tp1.lo[i]), as.character(tp1.c[i]), as.character(tp1.p[i]), as.character(tp1.s[i])),
    c(as.character(tp2.lo[i]), as.character(tp2.c[i]), as.character(tp2.p[i]), as.character(tp2.s[i])),
    c(as.character(tp3.lo[i]), as.character(tp3.c[i]), as.character(tp3.p[i]), as.character(tp3.s[i])),
    c(as.character(tp4.lo[i]), as.character(tp4.c[i]), as.character(tp4.p[i]), as.character(tp4.s[i])),
    c(as.character(tp5.lo[i]), as.character(tp5.c[i]), as.character(tp5.p[i]), as.character(tp5.s[i])),
    c(as.character(tp6.lo[i]), as.character(tp6.c[i]), as.character(tp6.p[i]), as.character(tp6.s[i]))
  )
  colnames(GMM[[i]])<-c("Indiv", "clus","FCvCaPval","status")
  rownames(GMM[[i]])<-c("tp_2min", "tp_10min", "tp_20min", "tp_30min", "tp_60min", "tp_90min")
}

GMM.wFC<-vector("list", dim(wtFC)[1])
names(GMM.wFC)<-rownames(wtFC)
for(i in 1:length(GMM.wFC)){
  GMM.wFC[[i]]<-rbind(
    c(as.character(tp1.lo[i]), as.character(tp1.c[i]), as.character(tp1.p[i]), as.character(tp1.s[i]), as.character(numData.Avg[i, 1])),
    c(as.character(tp2.lo[i]), as.character(tp2.c[i]), as.character(tp2.p[i]), as.character(tp2.s[i]), as.character(numData.Avg[i, 2])),
    c(as.character(tp3.lo[i]), as.character(tp3.c[i]), as.character(tp3.p[i]), as.character(tp3.s[i]), as.character(numData.Avg[i, 3])),
    c(as.character(tp4.lo[i]), as.character(tp4.c[i]), as.character(tp4.p[i]), as.character(tp4.s[i]), as.character(numData.Avg[i, 4])),
    c(as.character(tp5.lo[i]), as.character(tp5.c[i]), as.character(tp5.p[i]), as.character(tp5.s[i]), as.character(numData.Avg[i, 5])),
    c(as.character(tp6.lo[i]), as.character(tp6.c[i]), as.character(tp6.p[i]), as.character(tp6.s[i]), as.character(numData.Avg[i, 6]))
  )
  colnames(GMM.wFC[[i]])<-c("Indiv", "clus", "FCvCaPval", "status", "FCvC")
  rownames(GMM.wFC[[i]])<-c("tp_2min", "tp_10min", "tp_20min", "tp_30min", "tp_60min", "tp_90min")
}

###
# Saving GMM object as a list

GMM.ID <- data.IDmap
colnames(GMM.ID) <- c("dataID", "UPID", "S.cc")
GMM.ID <- as.data.frame(GMM.ID)

save(list=c("GMM.ID", "GMM","GMM.wFC"), file="dataGMM_A2058.RData")
