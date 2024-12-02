#### common DEGs ####

common_DEGs <- function(list_DEGs, names = c("sAAb, mAAb")){
  
  if (length(list_DEGs) == 2){
    
    DEGs_1 <- list_DEGs[[1]]
    DEGs_2 <- list_DEGs[[2]]
    
    overlap <- intersect(DEGs_1$ProbesetID, 
                         DEGs_2$ProbesetID)
    
    
    
    DEGs_1 <- DEGs_1 %>% mutate(diffexpressed = case_when(
      logFC > 0  ~ 'UP',
      logFC < 0  ~ 'DOWN'
    ))
    
    DEGs_1$Donor_type_expression <- rep(paste(names[1]," only", sep = ""), times = nrow(DEGs_1))
    
    
    
    DEGs_2 <- DEGs_2 %>% mutate(diffexpressed = case_when(
      logFC > 0  ~ 'UP',
      logFC < 0  ~ 'DOWN'
    ))
    
    DEGs_2$Donor_type_expression <- rep(paste(names[2]," only", sep = ""), times = nrow(DEGs_2))
    
    ## sAAb table
    for (i in 1:nrow(DEGs_1)){
      
      if(DEGs_1$ProbesetID[i] %in% overlap ){
        
        DEGs_1$Donor_type_expression[i] <- paste(names[1], " and ", names[2], sep = "")
      }
      
    }
    
    #maAb table
    for (i in 1:nrow(DEGs_2)){
      
      if(DEGs_2$ProbesetID[i] %in% overlap ){
        
        DEGs_2$Donor_type_expression[i] <- paste(names[1], " and ", names[2], sep = "")
      }
      
    }
    
    ## overlap table
    
    DEGs_1_common <- DEGs_1[,c(1,2,3,7,9,10,11)] %>% dplyr::filter(Donor_type_expression == paste(names[1], " and ", names[2], sep = ""))
    DEGs_1_common <- DEGs_1_common[order(DEGs_1_common[,1]),]
    colnames(DEGs_1_common) <- unlist(sapply(colnames(DEGs_1_common), paste, "_1", sep=""))
    
    DEGs_2_common <- DEGs_2[,c(1,2,3,7,9,10,11)] %>% dplyr::filter(Donor_type_expression == paste(names[1], " and ", names[2], sep = ""))
    DEGs_2_common <- DEGs_2_common[order(DEGs_2_common[,1]),]
    colnames(DEGs_2_common) <- unlist(sapply(colnames(DEGs_2_common), paste, "_2", sep=""))
    
    if(all(DEGs_1_common[,1] == DEGs_2_common[,1])){
      print("All rows matched")
      
      overlap_DEGs <- bind_cols( DEGs_1_common, DEGs_2_common )
      
    }else print("Rows do not match")
    
    overlap_DEGs <- overlap_DEGs %>% mutate(Direction = case_when(
      diffexpressed_1 == diffexpressed_2 ~ "Same dir",
      diffexpressed_1 != diffexpressed_2 ~ "Opposite"))
    
    overlap_DEGs <- overlap_DEGs[,c(1,5,2,3,10,15,14,4,11)]
    
    results<- list(DEGs_1 = DEGs_1,DEGs_2 = DEGs_2, Overlap = overlap_DEGs)
    
    
  } else{
    
    DEGs_1 <- list_DEGs[[1]]
    DEGs_2 <- list_DEGs[[2]]
    DEGs_3 <- list_DEGs[[3]]
    
    overlap <- intersect(intersect(DEGs_1$ProbesetID, 
                                   DEGs_2$ProbesetID),
                         DEGs_3$ProbesetID)
    
    
    
    DEGs_1 <- DEGs_1 %>% mutate(diffexpressed = case_when(
      logFC > 0  ~ 'UP',
      logFC < 0  ~ 'DOWN'
    ))
    
    DEGs_1$Donor_type_expression <- rep(paste(names[1]," only", sep = ""), times = nrow(DEGs_1))
    
    
    
    DEGs_2 <- DEGs_2 %>% mutate(diffexpressed = case_when(
      logFC > 0  ~ 'UP',
      logFC < 0  ~ 'DOWN'
    ))
    
    DEGs_2$Donor_type_expression <- rep(paste(names[2]," only", sep = ""), times = nrow(DEGs_2))
    
    DEGs_3 <- DEGs_3 %>% mutate(diffexpressed = case_when(
      logFC > 0  ~ 'UP',
      logFC < 0  ~ 'DOWN'
    ))
    
    DEGs_3$Donor_type_expression <- rep(paste(names[3]," only", sep = ""), times = nrow(DEGs_3))
    
    ## sAAb table
    for (i in 1:nrow(DEGs_1)){
      
      if(DEGs_1$ProbesetID[i] %in% overlap ){
        
        DEGs_1$Donor_type_expression[i] <- paste(names[1], ", ",  names[2], " and ", names[3], sep = "")
      }
      
    }
    
    #maAb table
    for (i in 1:nrow(DEGs_2)){
      
      if(DEGs_2$ProbesetID[i] %in% overlap ){
        
        DEGs_2$Donor_type_expression[i] <- paste(names[1], ", ",  names[2], " and ", names[3], sep = "")
      }
      
    }
    
    #3
    for (i in 1:nrow(DEGs_3)){
      
      if(DEGs_3$ProbesetID[i] %in% overlap ){
        
        DEGs_3$Donor_type_expression[i] <- paste(names[1], ", ",  names[2], " and ", names[3],  sep = "")
      }
      
    }
    
    ## overlap table
    
    DEGs_1_common <- DEGs_1[,c(1,2,3,7,9,10,11)] %>% dplyr::filter(Donor_type_expression == paste(names[1], ", ",  names[2], " and ", names[3], sep = ""))
    DEGs_1_common <- DEGs_1_common[order(DEGs_1_common[,1]),]
    colnames(DEGs_1_common) <- unlist(sapply(colnames(DEGs_1_common), paste, "_1", sep=""))
    
    DEGs_2_common <- DEGs_2[,c(1,2,3,7,9,10,11)] %>% dplyr::filter(Donor_type_expression == paste(names[1], ", ",  names[2], " and ", names[3], sep = ""))
    DEGs_2_common <- DEGs_2_common[order(DEGs_2_common[,1]),]
    colnames(DEGs_2_common) <- unlist(sapply(colnames(DEGs_2_common), paste, "_2", sep=""))
    
    DEGs_3_common <- DEGs_3[,c(1,2,3,7,9,10,11)] %>% dplyr::filter(Donor_type_expression == paste(names[1], ", ",  names[2], " and ", names[3], sep = ""))
    DEGs_3_common <- DEGs_3_common[order(DEGs_3_common[,1]),]
    colnames(DEGs_3_common) <- unlist(sapply(colnames(DEGs_3_common), paste, "_3", sep=""))
    
    
    if(all(DEGs_1_common[,1] == DEGs_2_common[,1]) && all(DEGs_2_common[,1] == DEGs_3_common[,1])){
      print("All rows matched")
      
      overlap_DEGs <- bind_cols( DEGs_1_common, DEGs_2_common, DEGs_3_common )
      
    }else print("Rows do not match")
    
    overlap_DEGs <- overlap_DEGs %>% mutate(Direction = case_when(
      (diffexpressed_1 == diffexpressed_2 & diffexpressed_2 == diffexpressed_3) ~ "Same dir",
      diffexpressed_1 != diffexpressed_2 ~ "Changing dir"))
    
    overlap_DEGs <- overlap_DEGs[,c(1,5,2,3,10, 17, 15,14,4,11, 18 )]
    
    
    results<- list(DEGs_1 = DEGs_1,DEGs_2 = DEGs_2, DEGs_3 = DEGs_3, Overlap = overlap_DEGs)
  }
  
  
  
  
  return(results)
  
}
