human_to_mouse_genes <- function(human_table ){
  
  ### retrieve ENTREZ IDs for human genes
  ## map ensembl gene ids
  
  print("ENSEMBL homolog search start")
  
  #results1 <-  getBM(c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), "hgnc_symbol", human_table$GeneName,mart=human, uniqueRows = TRUE)
  results <- getBM(c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), "ensembl_gene_id", human_table$ensembl_gene_ID,mart=human, uniqueRows = TRUE)
  #results <- dplyr::union(results1, results2)
  
  print("Entrez id mapping done")
  
  ## homolog search with ensembl
  
  print("Mouse honmolog search start")
  results_ensembl <- getBM(c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name"), "ensembl_gene_id", results$ensembl_gene_id,mart=human, uniqueRows = TRUE)
  
  ## getting human entrez
  entrezgene_id <- mapIds(org.Hs.eg.db, results_ensembl$ensembl_gene_id, "ENTREZID","ENSEMBL")
  entrezgene_id <- do.call(c, lapply(entrezgene_id, function(x) if(is.null(x)) return(NA) else return(x)))
  
  results_ensembl$entrezgene_id <- as.numeric(entrezgene_id)
  
  #adding hgnc symbol to results ensembl
  
  for(i in 1:nrow(results_ensembl)){
    
    for(j in 1:nrow(results)){
      
      if(results_ensembl$ensembl_gene_id[i] == results$ensembl_gene_id[j]){
        
        results_ensembl$hgnc_symbol[i] <- results$hgnc_symbol[j]
        
        
      }
      
    }
    
    
  }
  
  print("ENSEMBL homolog search done")
  
  ## mapping with NCBI
  print("NCBI Homolog mapping start")
  mapped <-select(Orthology.eg.db, results_ensembl$entrezgene_id, "Mus.musculus","Homo.sapiens")
  print("NCBI Homolog mapping done")
  
  # getting symbol and ensembl for mouse
  mouse_mapped <- mapIds(org.Mm.eg.db, as.character(mapped[,2]), "ENSEMBL","ENTREZID")
  mouse_mapped2 <- mapIds(org.Mm.eg.db, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  
  #unlisting previous results
  mouse_mapped <- do.call(c, lapply(mouse_mapped, function(x) if(is.null(x)) return(NA) else return(x)))
  mouse_mapped2 <- do.call(c, lapply(mouse_mapped2, function(x) if(is.null(x)) return(NA) else return(x)))
  
  # human ensembl
  # human_mapped <- mapIds(org.Hs.eg.db, as.character(mapped[,1]), "ENSEMBL","ENTREZID")
  # human_mapped2 <- mapIds(org.Hs.eg.db, as.character(mapped[,1]), "SYMBOL","ENTREZID")
  
  # human_mapped<- do.call(c, lapply(human_mapped, function(x) if(is.null(x)) return(NA) else return(x)))
  # human_mapped2<- do.call(c, lapply(human_mapped2, function(x) if(is.null(x)) return(NA) else return(x)))
  
  #binding all in one table
  final_mapped <- cbind(results_ensembl$ensembl_gene_id, results_ensembl$hgnc_symbol,results_ensembl$entrezgene_id,mapped[,2], mouse_mapped2, mouse_mapped)
  
  colnames(final_mapped) <- c("ensembl_gene_id","hgnc_symbol","entrezgene_id", "Mm_entrezgene_id","mmusculus_homolog_associated_gene_name", "mmusculus_homolog_ensembl_gene")
  
  final_mapped <- as.data.frame(final_mapped)
  final_mapped$entrezgene_id <- as.numeric(final_mapped$entrezgene_id)
  final_mapped$Mm_entrezgene_id <- as.numeric(final_mapped$Mm_entrezgene_id)
  
  
  final_mapped_Ms <- full_join(results_ensembl, final_mapped)
  final_mapped_Ms <- unique(final_mapped_Ms)
  
  print("Homolog mapping done")
  
  ## sorting DEGs table with mapped table
  print("Sorting mapped table with DEG table")
  final_mapped_sorted <- final_mapped_Ms[order(match(final_mapped_Ms[,"ensembl_gene_id"], human_table$ensembl_gene_ID)) ,]
  
  
  
  ## adding CD3 stain label 
  print("CD3 label start")
  
  label_gene <- NULL
  for(i in 1:nrow(final_mapped_sorted)){
    
    for (j in 1:nrow(human_table)){
      
      if(grepl("^$", final_mapped_sorted[,"hgnc_symbol"][i]) || is.na(final_mapped_sorted[,"hgnc_symbol"][i])){
        
        if (final_mapped_sorted[,"ensembl_gene_id"][i] == human_table$ensembl_gene_ID[j]){
          label_eval <- human_table$Label[j]
          print(i)
        }
        
      } else if (final_mapped_sorted[,"hgnc_symbol"][i] == human_table$GeneName[j]){
        label_eval <- human_table$Label[j]
        print(i)
      }
      
      
    }
    label_gene <- c(label_gene, label_eval)
  }
  
  final_mapped_sorted$Label<- label_gene
  
  return(final_mapped_sorted)
  
}