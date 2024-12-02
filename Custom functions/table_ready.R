## average and sd gene expression subsets for plotting curves

table_ready <- function (gene_expression, genes_of_interest, groups = c("ND", "sAAb", "mAAb", "T1D"), expression_values = F){
  
  mean_sd_expression <- NULL
  
  if (nrow(genes_of_interest) == 1){
    genes_interest_expression <- t(as.matrix(gene_expression[rownames(gene_expression) %in% rownames(genes_of_interest),]))
    rownames(genes_interest_expression) <- rownames(genes_of_interest)
    
    for(i in 1:length(groups)){
      a<- NULL
      a <- genes_interest_expression[,grepl(groups[i],colnames(genes_interest_expression))]
      probe_id <- rownames(genes_interest_expression)
      average_expression <- mean(a) 
      sd_gene_expression <- sd(a)
      
      serror_gene_expression <- std.error(a)
      
      average_sd_expression <- data.frame(probe_id = probe_id, genesymbol=genes_of_interest$genesymbol, mRNA_accession = genes_of_interest$mRNA_Accession, 
                                          Clinical_phenotype = rep(groups[i], times=length(rownames(genes_interest_expression))), Expression = average_expression, Sd = sd_gene_expression,
                                          Std_error = serror_gene_expression)
      mean_sd_expression <- rbind(mean_sd_expression, average_sd_expression)
      
    }
    
  }else{
    
    genes_interest_expression <- gene_expression[rownames(gene_expression) %in% rownames(genes_of_interest),]
    genes_interest_expression <- genes_interest_expression[order(match(rownames(genes_interest_expression), rownames(genes_of_interest))),]
    
    for(i in 1:length(groups)){
      a<- NULL
      a <- genes_interest_expression[,grepl(groups[i],colnames(genes_interest_expression))]
      probe_id <- rownames(a)
      average_expression <- apply(a, 1, mean) 
      sd_gene_expression <- apply(a,1, sd)
      
      serror_gene_expression <- apply(a, 1, std.error)
      
      average_sd_expression <- data.frame(probe_id = probe_id, genesymbol=genes_of_interest$genesymbol, mRNA_accession = genes_of_interest$mRNA_Accession, 
                                          Clinical_phenotype = rep(groups[i], times=length(rownames(a))), Expression = average_expression, Sd = sd_gene_expression,
                                          Std_error = serror_gene_expression)
      mean_sd_expression <- rbind(mean_sd_expression, average_sd_expression)
      
    }
  }
  
  
  
  
  #colnames(mean_expression) <- groups
  rownames(mean_sd_expression) <- NULL
  mean_sd_expression$Clinical_phenotype <- factor(mean_sd_expression$Clinical_phenotype, levels=unique(mean_sd_expression$Clinical_phenotype))
  
  if (expression_values == T){
    genes_interest_expression <- cbind(genes_interest_expression, genes_of_interest$genesymbol)
    results <- list(genes_interest_expression, mean_sd_expression)
    return(results)
  }else return(mean_sd_expression)
  
  
}