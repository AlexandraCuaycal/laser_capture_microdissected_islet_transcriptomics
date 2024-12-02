matching_mouse_human <- function(mouse_DEGs, human_to_mouse_table, human_table){
  
  mouse_DEGs_id <- mouse_DEGs$ensembl_gene_id
  
  matched_DEGs_human <- human_to_mouse_table[human_to_mouse_table$mmusculus_homolog_ensembl_gene %in% mouse_DEGs_id, ]
  
  matched_DEGs_human_table <- human_table[human_table$ensembl_gene_ID %in% matched_DEGs_human$ensembl_gene_id, ]
  
  matched_DEGs_mice <- mouse_DEGs[mouse_DEGs$ensembl_gene_id %in% matched_DEGs_human$mmusculus_homolog_ensembl_gene, ]
  
  label_table <- matched_DEGs_human[,c(1,2,3,5,7)]
  colnames(label_table)[c(1,2)]<-c( "Hs_ensembl_gene_id","ensembl_gene_id")
  
  label_table <- label_table[!duplicated(label_table[c("ensembl_gene_id")]),]
  label_table <- label_table[label_table$ensembl_gene_id %in% mouse_DEGs_id, ]
  
  matched_DEGs_mice <-dplyr::left_join(matched_DEGs_mice, label_table)
  
  results <- list(Matched_in_human=matched_DEGs_human_table, Matched_in_mice = matched_DEGs_mice, Matched_table=matched_DEGs_human)
  
  return(results)
  
}