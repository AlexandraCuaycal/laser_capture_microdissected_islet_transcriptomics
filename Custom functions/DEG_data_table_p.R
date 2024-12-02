### Custom functions to save / plot results

## function to get protein expression data with adj p val = 0.05

DEP_data_table <- function(DEP_object, adj.p.val = 0.05, logFCcutoff = 0.5, folder=folder, contrast_to_get = contrast_to_get){
  
  DEPs_by_contrast<- list()
  
  
  for(i in 1:length(contrast_to_get)){
    
    index_sig_DEPs <- which(DEP_object[[i]]$adj.P.Val < adj.p.val)
    index_logFC <- which(abs(DEP_object[[i]]$logFC) > logFCcutoff )
    index_DEPs <- intersect(index_sig_DEPs, index_logFC)
    sig_DEPs_table <- DEP_object[[i]][index_DEPs,]

    DEPs_by_contrast[[i]] <- sig_DEPs_table[sort(sig_DEPs_table$adj.P.Val, index.return=TRUE)$ix,]
    DEPs_by_contrast[[i]] <- DEPs_by_contrast[[i]][order(DEPs_by_contrast[[i]]$logFC, decreasing = TRUE),]
    
    if (!is.null(folder)) {
      
      dir.create(folder)
      
      write.csv(DEPs_by_contrast[[i]], file = paste(folder, "/", contrast_to_get[i], ".csv", sep=""))
    }
    
  }
  
  names(DEPs_by_contrast) <- contrast_to_get
  return(DEPs_by_contrast)
  
}