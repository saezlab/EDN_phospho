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
# 13:00 23/03/2018
# This function is used to build a matrix containing the scored assigned for
# each of the measured nodes

buildDataMatrix <- function(dataGMM = dataGMM, pknList = pknList, targets, experiments){
  
  tg <- unlist(targets)
  
  allSpecies <- pknList@species
  dataMatrix <- matrix(, nrow = length(unlist(experiments)), ncol = length(allSpecies))
  
  #colnames
  cNames <- c()
  kk <- which(allSpecies%in%tg)
  for(i in 1:length(kk)){
    cNames <- c(cNames, paste("TG:", allSpecies[kk[i]], sep = ""))
  }
  
  idxDN <- which(!(allSpecies%in%dataGMM@IDmap$S.cc))
  idxDN <- setdiff(idxDN, kk)
  for(i in 1:length(idxDN)){
    cNames <- c(cNames, paste("DN:", allSpecies[idxDN[i]], sep = ""))
  }
  
  dataMatrix[, 1:length(cNames)] <- as.numeric(0)
  
  idxDS <- which(allSpecies%in%dataGMM@IDmap$S.cc)
  ds <- length(cNames)
  for(i in 1:length(idxDS)){
    cNames <- c(cNames, paste("DS:", allSpecies[idxDS[i]], sep = ""))
  }
  
  colnames(dataMatrix) <- cNames
  
  #set values on the matrix
  for(i in 1:length(idxDS)){
    cpSite <- allSpecies[idxDS[i]]
    dataSite <- as.character(dataGMM@IDmap$dataID[which(dataGMM@IDmap$S.cc==cpSite)])
    for(j in 1:length(dataSite)){
      dataID <- which(names(dataGMM@res)==dataSite[j])
      if(length(dataID) > 0){
        for(dd in 1:length(dataID)){
          cnt <- 1
          for(k in 1:length(experiments)){
            for(l in 1:length(experiments[[k]])){
              if(!is.na(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 4])){
                if(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 4] == "OK"){
                  dataMatrix[cnt, ds+i] <- as.numeric(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 1])
                  cnt <- cnt + 1
                }
              }
            }
          }
        }
      }
    }
  }
  dataMatrix[is.na(dataMatrix)] <- as.numeric(0)
  
  #indeces
  tgID <- c()
  dnID <- c()
  dsID <- c()
  for(i in 1:ncol(dataMatrix)){
    
    if(strsplit(colnames(dataMatrix)[i], split = ":")[[1]][1]=="TG"){
      tgID <- c(tgID, i)
    }
    
    if(sum(dataMatrix[, i])==0){
      dnID <- c(dnID, i)
    }
    dnID <- dnID[!(dnID%in%tgID)]
    
    if(sum(dataMatrix[, i]!=0)){
      dsID <- c(dsID, i)
    }
    
  }
  
  #species
  species <- c()
  for(i in 1:ncol(dataMatrix)){
    
    species <- c(species, strsplit(colnames(dataMatrix)[i], split = ":")[[1]][2])
    
  }
  
  toRem <- c()
  kk <- which(colSums(dataMatrix)==0)
  for(i in 1:length(kk)){
    
    if(strsplit(names(kk[i]), split = ":")[[1]][1] == "DS"){
      
      toRem <- c(toRem, kk[i])
      
    }
    
  }
  
  if(length(toRem) > 0){
    
    dataMatrix[, toRem] <- as.numeric(1)
    
  }
  
  # dataMatrix <- as.matrix(t(dataMatrix))
  
  res <- list(dataMatrix=dataMatrix, tgID=tgID, dnID=dnID, dsID=dsID, species=species)
  
  return(res)
}