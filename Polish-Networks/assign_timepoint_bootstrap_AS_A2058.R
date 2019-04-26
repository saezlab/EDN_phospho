load(file = "../Analysis/resList_A2058.RData")

RESULT <- matrix(data = , nrow = 1, ncol = 3)
colnames(RESULT) <- c("Source", "Entry", "Target")

for(ii in 1:length(resList)){
  
  resultSIF_tp1 <- as.data.frame(resList[[ii]][[1]])
  resultSIF_tp1$Source <- as.character(resultSIF_tp1$Source)
  resultSIF_tp1$Interaction <- as.character(resultSIF_tp1$Interaction)
  resultSIF_tp1$Target <- as.character(resultSIF_tp1$Target)
  resultSIF_tp1$X <-paste(resultSIF_tp1[,1], resultSIF_tp1[,3], sep = "|")
  
  resultSIF_tp2 <- as.data.frame(resList[[ii]][[2]])
  resultSIF_tp2$Source <- as.character(resultSIF_tp2$Source)
  resultSIF_tp2$Interaction <- as.character(resultSIF_tp2$Interaction)
  resultSIF_tp2$Target <- as.character(resultSIF_tp2$Target)
  resultSIF_tp2$X <-paste(resultSIF_tp2[,1], resultSIF_tp2[,3], sep = "|")
  
  resultSIF_tp3 <- as.data.frame(resList[[ii]][[3]])
  resultSIF_tp3$Source <- as.character(resultSIF_tp3$Source)
  resultSIF_tp3$Interaction <- as.character(resultSIF_tp3$Interaction)
  resultSIF_tp3$Target <- as.character(resultSIF_tp3$Target)
  resultSIF_tp3$X <-paste(resultSIF_tp3[,1], resultSIF_tp3[,3], sep = "|")
  
  resultSIF_tp4 <- as.data.frame(resList[[ii]][[4]])
  resultSIF_tp4$Source <- as.character(resultSIF_tp4$Source)
  resultSIF_tp4$Interaction <- as.character(resultSIF_tp4$Interaction)
  resultSIF_tp4$Target <- as.character(resultSIF_tp4$Target)
  resultSIF_tp4$X <-paste(resultSIF_tp4[,1], resultSIF_tp4[,3], sep = "|")
  
  resultSIF_tp5 <- as.data.frame(resList[[ii]][[5]])
  resultSIF_tp5$Source <- as.character(resultSIF_tp5$Source)
  resultSIF_tp5$Interaction <- as.character(resultSIF_tp5$Interaction)
  resultSIF_tp5$Target <- as.character(resultSIF_tp5$Target)
  resultSIF_tp5$X <-paste(resultSIF_tp5[,1], resultSIF_tp5[,3], sep = "|")
  
  resultSIF_tp6 <- as.data.frame(resList[[ii]][[6]])
  resultSIF_tp6$Source <- as.character(resultSIF_tp6$Source)
  resultSIF_tp6$Interaction <- as.character(resultSIF_tp6$Interaction)
  resultSIF_tp6$Target <- as.character(resultSIF_tp6$Target)
  resultSIF_tp6$X <-paste(resultSIF_tp6[,1], resultSIF_tp6[,3], sep = "|")
  
  resultSIF_tp6$Y<-"90min"
  int<-intersect(resultSIF_tp6$X, resultSIF_tp5$X)
  resultSIF_tp6$Y[which(match(resultSIF_tp6$X, int) >0)] <- "60min"
  
  int<-intersect(resultSIF_tp5$X, resultSIF_tp4$X)
  resultSIF_tp6$Y[which(match(resultSIF_tp6$X, int) >0)] <- "30min"
  
  int<-intersect(resultSIF_tp4$X, resultSIF_tp3$X)
  resultSIF_tp6$Y[which(match(resultSIF_tp6$X, int) >0)] <- "20min"
  
  int<-intersect(resultSIF_tp3$X, resultSIF_tp2$X)
  resultSIF_tp6$Y[which(match(resultSIF_tp6$X, int) >0)] <- "10min"
  
  int<-intersect(resultSIF_tp2$X, resultSIF_tp1$X)
  resultSIF_tp6$Y[which(match(resultSIF_tp6$X, int) >0)] <- "2min"
  
  result<-cbind(resultSIF_tp6[,1],resultSIF_tp6[,5],resultSIF_tp6[,3])
  colnames(result)<-c("Source","Entry","Target")
  
  RESULT <- unique(rbind(RESULT, result))
  
}

RESULT <- RESULT[-1, ]

attribSIF <- matrix(data = , nrow = nrow(RESULT), ncol = 4)
colnames(attribSIF) <- c(colnames(RESULT), "weight")

attribSIF[, 1:3] <- RESULT

for(i in 1:nrow(attribSIF)){
  
  cnt <- 0
  
  #
  if(attribSIF[i, 2]=="2min"){
    
    for(j in 1:length(resList)){
      
      idx1 <- which(resList[[j]][[1]][, 1]==attribSIF[i, 1])
      idx2 <- which(resList[[j]][[1]][, 3]==attribSIF[i, 3])
      
      if((length(idx1) > 0) && (length(idx2) > 0)){
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          cnt <- cnt + 1
          
        }
        
      }
      
    }
    
  }
  
  #
  if(attribSIF[i, 2]=="10min"){
    
    for(j in 1:length(resList)){
      
      idx1 <- which(resList[[j]][[2]][, 1]==attribSIF[i, 1])
      idx2 <- which(resList[[j]][[2]][, 3]==attribSIF[i, 3])
      
      if((length(idx1) > 0) && (length(idx2) > 0)){
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          cnt <- cnt + 1
          
        }
        
      }
      
    }
    
  }
  
  #
  if(attribSIF[i, 2]=="20min"){
    
    for(j in 1:length(resList)){
      
      idx1 <- which(resList[[j]][[3]][, 1]==attribSIF[i, 1])
      idx2 <- which(resList[[j]][[3]][, 3]==attribSIF[i, 3])
      
      if((length(idx1) > 0) && (length(idx2) > 0)){
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          cnt <- cnt + 1
          
        }
        
      }
      
    }
    
  }
  
  
  #
  if(attribSIF[i, 2]=="30min"){
    
    for(j in 1:length(resList)){
      
      idx1 <- which(resList[[j]][[4]][, 1]==attribSIF[i, 1])
      idx2 <- which(resList[[j]][[4]][, 3]==attribSIF[i, 3])
      
      if((length(idx1) > 0) && (length(idx2) > 0)){
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          cnt <- cnt + 1
          
        }
        
      }
      
    }
    
  }
  
  #
  if(attribSIF[i, 2]=="60min"){
    
    for(j in 1:length(resList)){
      
      idx1 <- which(resList[[j]][[5]][, 1]==attribSIF[i, 1])
      idx2 <- which(resList[[j]][[5]][, 3]==attribSIF[i, 3])
      
      if((length(idx1) > 0) && (length(idx2) > 0)){
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          cnt <- cnt + 1
          
        }
        
      }
      
    }
    
  }
  
  #
  if(attribSIF[i, 2]=="90min"){
    
    for(j in 1:length(resList)){
      
      idx1 <- which(resList[[j]][[6]][, 1]==attribSIF[i, 1])
      idx2 <- which(resList[[j]][[6]][, 3]==attribSIF[i, 3])
      
      if((length(idx1) > 0) && (length(idx2) > 0)){
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          cnt <- cnt + 1
          
        }
        
      }
      
    }
    
  }
  
  attribSIF[i, 4] <- as.character(cnt)
  
}

##
uniqueAttribSIF <- matrix(data = , nrow = nrow(unique(attribSIF[, c(1, 3)])), ncol = 4)
colnames(uniqueAttribSIF) <- colnames(attribSIF)

uniqueAttribSIF[, c(1, 3)] <- unique(attribSIF[, c(1, 3)])
for(i in 1:nrow(uniqueAttribSIF)){
  
  idx1 <- which(attribSIF[, 1]==uniqueAttribSIF[i, 1])
  idx2 <- which(attribSIF[, 3]==uniqueAttribSIF[i, 3])
  
  idx <- intersect(x = idx1, y = idx2)
  
  if(length(idx)==1){
    
    uniqueAttribSIF[i, 2] <- attribSIF[idx, 2]
    uniqueAttribSIF[i, 4] <- attribSIF[idx, 4]
    
  } else {
    
    cntMaxIdx <- which(as.numeric(attribSIF[idx, 4])==max(as.numeric(attribSIF[idx, 4])))
    uniqueAttribSIF[i, 4] <- attribSIF[idx[cntMaxIdx], 4]
    #uniqueAttribSIF[i, 4] <- as.character(mean(as.numeric(attribSIF[idx, 4])))
     tp1<-attribSIF[idx[which(attribSIF[idx,2]=="2min")],]
     tp2<-attribSIF[idx[which(attribSIF[idx,2]=="10min")],]
     tp3<-attribSIF[idx[which(attribSIF[idx,2]=="20min")],]
     tp4<-attribSIF[idx[which(attribSIF[idx,2]=="30min")],]
     tp5<-attribSIF[idx[which(attribSIF[idx,2]=="60min")],]
     tp6<-attribSIF[idx[which(attribSIF[idx,2]=="90min")],]
     if(length(tp1)<1){tp1<-c(1,1,1,1)}
     if(length(tp2)<1){tp2<-c(1,1,1,1)}
     if(length(tp3)<1){tp3<-c(1,1,1,1)}
     if(length(tp4)<1){tp4<-c(1,1,1,1)}
     if(length(tp5)<1){tp5<-c(1,1,1,1)}
     if(length(tp6)<1){tp6<-c(1,1,1,1)}
     
     if(as.numeric(tp1[4]) > 20){uniqueAttribSIF[i, 2]<-tp1[2]} 
     else {
      if((as.numeric(tp2[4]) > 20)){uniqueAttribSIF[i, 2]<-tp2[2]}
      else {
        if((as.numeric(tp3[4]) > 20)){uniqueAttribSIF[i, 2]<-tp3[2]}
      else { 
        if((as.numeric(tp4[4]) > 20)){uniqueAttribSIF[i, 2]<-tp4[2]}
        else { 
          if((as.numeric(tp5[4]) > 20)){uniqueAttribSIF[i, 2]<-tp5[2]}
        else {uniqueAttribSIF[i, 2]<-tp6[2]}
        
       }
        }
      }
     }
  }
  
}

write.table(x = uniqueAttribSIF, file = "inst/tpBootstrapAnalysis_A2058.txt", quote = FALSE, sep = "\t", row.names = FALSE)
