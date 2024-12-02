### Custom functions to save / plot results

## function to get gene expression data with adj p val = 0.001

DEG_data_table <- function(DEG_object=DEG_object, adj.p.val = 0.001, logFCcutoff = 0.5, folder=folder, contrast_to_get = contrast_to_get){
  
  DEGs_by_contrast<- list()
  for(i in 1:length(contrast_to_get)){
    
    index_sig_DEGs <- which(DEG_object$resTable[[i]]$adj.P.Val < adj.p.val)
    index_logFC <- which(abs(DEG_object$resTable[[i]]$logFC) > logFCcutoff )
    index_DEGs <- intersect(index_sig_DEGs, index_logFC)
    sig_DEGs_table <- DEG_object$resTable[[i]][index_DEGs,]
    ensembl_IDs <- ensembl_gene_id[rownames(ensembl_gene_id) %in% sig_DEGs_table$ProbesetID , ]
    ensembl_IDs <-ensembl_IDs[ match(rownames(ensembl_IDs),sig_DEGs_table$ProbesetID ),]
    sig_DEGs_table$ensembl_gene_ID <- ensembl_IDs$ensembl_gene_id
    DEGs_by_contrast[[i]] <- sig_DEGs_table[sort(sig_DEGs_table$adj.P.Val, index.return=TRUE)$ix,]
    DEGs_by_contrast[[i]] <- DEGs_by_contrast[[i]][order(DEGs_by_contrast[[i]]$logFC, decreasing = TRUE),]
    
    if (!is.null(folder)) write.csv(DEGs_by_contrast[[i]], file = paste(folder, "/", contrast_to_get[i], ".csv", sep=""))
    
  }
  
  names(DEGs_by_contrast) <- contrast_to_get
  return(DEGs_by_contrast)
  
}