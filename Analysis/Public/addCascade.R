addCascade <- function(cascade=cascade, allD=allD, pknList=pknList){
  
  interactions <- c()
  for(ii in 1:(length(cascade)-1)){
    
    for(i in 1:length(cascade[[ii]])){
      
      kin <- cascade[[ii]][i]
      
      for(j in 1:length(cascade[[ii+1]])){
        
        sbstr <- cascade[[ii+1]][j]
        
        # idx <- intersect(which(allD$K.ID==kin))
        interactions <- c(interactions, intersect(which(allD$K.ID==kin), which(allD$S.ID==sbstr)))
        
      }
      
    }
    
  }
  
  toBind <- matrix(data = , nrow = length(interactions), ncol = ncol(pknList@interactions))
  colnames(toBind) <- colnames(pknList@interactions)
  for(i in 1:length(interactions)){
    
    toBind[i, 1] <- allD[interactions[i], 1]
    toBind[i, 2] <- allD[interactions[i], 2]
    toBind[i, 3] <- allD[interactions[i], 3]
    toBind[i, 4] <- allD[interactions[i], 4]
    toBind[i, 5] <- allD[interactions[i], 5]
    toBind[i, 6] <- allD[interactions[i], 6]
    toBind[i, 7] <- allD[interactions[i], 7]
    toBind[i, 8] <- allD[interactions[i], 8]
    toBind[i, 9] <- "1"
    
  }
  
  toBind <- as.data.frame(toBind)
  
  toBind$S.AC <- as.character(toBind$S.AC)
  toBind$S.ID <- as.character(toBind$S.ID)
  toBind$K.AC <- as.character(toBind$K.AC)
  toBind$K.ID <- as.character(toBind$K.ID)
  toBind$res <- as.character(toBind$res)
  toBind$pos <- as.character(toBind$pos)
  toBind$SID <- as.character(toBind$SID)
  toBind$S.cc <- as.character(toBind$S.cc)
  toBind$ntag <- as.character(toBind$ntag)
  
  pknList@interactions <- rbind(toBind, pknList@interactions)
  
  bb <- toBind
  
  for(i in 1:nrow(bb)){
    
    toBind <- matrix(data = , nrow = 1, ncol = ncol(pknList@interactions))
    colnames(toBind) <- colnames(pknList@interactions)
    
    idx <- which(is.na(pknList@interactions[1]))
    
    if((pknList@interactions$S.cc[i]%in%pknList@interactions$K.ID[idx])==FALSE){
      
      toBind[1, 4] <- as.character(pknList@interactions$S.cc[i])
      toBind[1, 7] <- paste0("i", as.numeric(substr(x = pknList@interactions$SID[nrow(pknList@interactions)], start = 2, stop = nchar(pknList@interactions$SID[nrow(pknList@interactions)])))+1)
      spl <- strsplit(x = as.character(pknList@interactions$S.cc[i]), split = "_")[[1]]
      toBind[1, 8] <- spl[1]
      for(j in 2:(length(spl)-1)){
        toBind[1, 8] <- paste0(toBind[1, 8], "_", spl[j])
      }
      toBind[1, 9] <- "1"
      
      toBind <- as.data.frame(toBind)
      
      toBind$K.ID <- as.character(toBind$K.ID)
      toBind$SID <- as.character(toBind$SID)
      toBind$S.cc <- as.character(toBind$S.cc)
      toBind$ntag <- as.character(toBind$ntag)
      
      pknList@interactions <- rbind(pknList@interactions, toBind)
      
    }
    
  }
  
  pkn <- pknList@interactions
  pkn <- pkn[, c(1, 2, 4, 5, 6, 8)]
  pkn <- unique(pkn)
  
  pknExp <- matrix(data = NA, nrow = nrow(pkn), ncol = 9)
  colnames(pknExp) <- colnames(pknList@interactions)
  
  cntE <- 1
  cntI <- 1
  
  for(i in 1:nrow(pknExp)){
    
    pknExp[i, 1] <- as.character(pkn$S.AC[i])
    pknExp[i, 2] <- as.character(pkn$S.ID[i])
    pknExp[i, 4] <- as.character(pkn$K.ID[i])
    pknExp[i, 5] <- as.character(pkn$res[i])
    pknExp[i, 6] <- as.character(pkn$pos[i])
    if(is.na(pknExp[i, 1])){
      pknExp[i, 7] <- paste0("i", cntI)
      cntI <- cntI+1
    } else {
      pknExp[i, 7] <- paste0("e", cntE)
      cntE <- cntE+1
    }
    pknExp[i, 8] <- as.character(pkn$S.cc[i])
    pknExp[i, 9] <- "1"
    pknExp[i, 3] <- as.character(pknList@interactions$K.AC)[which(pknList@interactions$K.ID==pknExp[i, 4])[1]]
    
  }
  
  pknExp <- as.data.frame(pknExp)
  
  pknExp$S.AC <- as.character(pknExp$S.AC)
  pknExp$S.ID <- as.character(pknExp$S.ID)
  pknExp$K.AC <- as.character(pknExp$K.AC)
  pknExp$K.ID <- as.character(pknExp$K.ID)
  pknExp$res <- as.character(pknExp$res)
  pknExp$pos <- as.character(pknExp$pos)
  pknExp$SID <- as.character(pknExp$SID)
  pknExp$S.cc <- as.character(pknExp$S.cc)
  pknExp$ntag <- as.character(pknExp$ntag)
  
  pknList@interactions <- pknExp
  pknList@species <- unique(c(pknList@interactions$K.ID, pknList@interactions$S.cc))
  
  return(pknList)
  
}