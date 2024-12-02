#############

# FUNCTION TO GET GENES FROM ENRICHMENTS ##########

genes_in_enrich <- function(enrichment, DEG_table, description, gene_description){
  
  ## getting column with genes
  genes_core_enrich <- enrichment@result %>% 
    dplyr::filter(Description %in% description)
  
  if(sum(colnames(genes_core_enrich) == "core_enrichment") == 1){
    
    genes_core_enrich <- dplyr::pull(genes_core_enrich, core_enrichment)
    
    ## separating if multiple columns with genes
    genes_enrichment <- NULL
    for( i in 1:length(genes_core_enrich)){
      
      genes <- unlist(strsplit(genes_core_enrich[i], "/", fixed = T))
      genes_enrichment <- c(genes_enrichment, genes)
      
    }
    genes_enrichment <- unique( genes_enrichment)
    
    ## getting probe id
    gene_description_enrich <- gene_description[ gene_description$ensembl_gene_id %in% genes_enrichment,]
    
  }else{
    genes_core_enrich <- dplyr::pull(genes_core_enrich, geneID)
    
    ## separating if multiple columns with genes
    genes_enrichment <- NULL
    for( i in 1:length(genes_core_enrich)){
      
      genes <- unlist(strsplit(genes_core_enrich[i], "/", fixed = T))
      genes_enrichment <- c(genes_enrichment, genes)
      
    }
    genes_enrichment <- unique( genes_enrichment)
    
    ## getting probe id
    gene_description_enrich <- gene_description[ gene_description$genesymbol %in% genes_enrichment,]
    
  }
  
  
  
  ## result with probe_id and genesymbol only
  genes_enrich_table <- gene_description_enrich %>% dplyr::select(genesymbol)
  genes_enrich_table$probe_id <- rownames(genes_enrich_table)
  
  #filtering with probe ids in the DEG table
  DEG_table_enrich <- DEG_table[DEG_table$ProbesetID %in% genes_enrich_table$probe_id,]
  
  colnames(DEG_table_enrich)[1:2] <- c("probe_id", "genesymbol")
  
  genes_unique <- unique(DEG_table_enrich$ensembl_gene_ID)
  
  index_max <- NULL
  
  for(i in 1:length(genes_unique)){
    
    table_subset <- DEG_table_enrich[DEG_table_enrich$ensembl_gene_ID == genes_unique[i],]
    max_val <- max( abs(table_subset$logFC) )
    
    if(all(abs(table_subset$logFC) ==  max_val)){
      index <-which(abs(table_subset$logFC) == max_val)[1]
      probe <- table_subset[index,1]
      index_f <- which(DEG_table_enrich$probe_id == probe)
    }else{
      index <- which(abs(table_subset$logFC) == max_val)[1]
      probe <- table_subset[index,1]
      index_f <- which(DEG_table_enrich$probe_id %in% probe)
    }
    
    index_max <- c(index_max, index_f)
    
  }
  
  DEG_table_enrich <- DEG_table_enrich[index_max,]
    
  
  return(DEG_table_enrich)
  
}