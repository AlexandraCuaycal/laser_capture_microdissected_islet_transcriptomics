## anovas and tukeys

tukey_tests <- function(gene_expression, genes_of_interest){
  
  phenotype <-NULL
  anovas<- list()
  tukeys <- list()
  labels<- NULL
  probes<- NULL
  probe_id<-NULL
  tukey_test_list <- NULL
  anova_dataset <- list()
  genes<- NULL
  comparisons <- NULL
  
  
  if (nrow(genes_of_interest) == 1){
    genes_interest_expression <- t(as.matrix(gene_expression[rownames(gene_expression) %in% rownames(genes_of_interest),]))
    rownames(genes_interest_expression) <- rownames(genes_of_interest)
    
    
  }else{
    
    genes_interest_expression <- gene_expression[rownames(gene_expression) %in% rownames(genes_of_interest),]
    genes_interest_expression <- genes_interest_expression[order(match(rownames(genes_interest_expression), rownames(genes_of_interest))),]
   
  }
  
  
  
  Clinical_phenotype <- strsplit(x = colnames(genes_interest_expression), split = "_")
  
  for(a in 1:length(Clinical_phenotype)){
    
    phenotype <- c(phenotype, Clinical_phenotype[[a]][1])
  }
  
  
  for (j in 1:length(rownames(genes_interest_expression))){
    
    
    anova_dataset[[j]] <- data.frame(probe_id = rep(rownames(genes_interest_expression)[j], times= length(colnames(genes_interest_expression))), 
                                     genesymbol=rep(genes_of_interest$genesymbol[j], times=length(colnames(genes_interest_expression))),
                                     Clinical_phenotype=phenotype, Expression = as.numeric(genes_interest_expression[j,]))
    
  }
  names(anova_dataset) <- rownames(genes_interest_expression)
  
  
  for (a in 1:length(anova_dataset)){
    
    anovas[[a]]<- aov(Expression ~ Clinical_phenotype, data = anova_dataset[[a]])
    tukeys[[a]] <- TukeyHSD(anovas[[a]])
    
    
    
    ## creating columns for tukey test tables
    
    labels <- c(labels, anova_dataset[[a]]$genesymbol[a])
    probes <- c(probes, anova_dataset[[a]]$probe_id[1])
    
    genes_n <- rep(labels[a], times = length(tukeys[[a]][[1]][,1]))
    genes<- c(genes, genes_n)
    
    probes_n <- rep(probes[a], times = length(tukeys[[a]][[1]][,1]))
    probe_id <- c(probe_id, probes_n)
    
    comp <- rownames(tukeys[[a]][[1]])
    comparisons <- c(comparisons, comp)
    
    tukey_test_list <- rbind(tukey_test_list, tukeys[[a]][[1]])
    
  }
  
  names(anovas) <- labels
  names(tukeys) <- labels
  tukey_test_list <- as.data.frame(tukey_test_list)
  
  
  
  ## tukey tests table
  
  tukey_list <- data.frame(probe_id = probe_id, genesymbol=genes, Comparison = comparisons, diff=tukey_test_list[,1], lwr = tukey_test_list[,2],
                           upr = tukey_test_list[,3], padj = tukey_test_list[,4]
  )
  
  
  tukey_list$padj_significance <- unlist(lapply(tukey_list$padj, function(x) {
    
    if(x >= 0.1) {
      x <- ""} else if (x > 0.05 && x < 0.1){
        x <- ""} else if(x <= 0.05 && x > 0.01){
          x <- "*"} else if (x <= 0.01 && x>0.001){
            x <- "**"}else if ( x <= 0.001 && x>0.0001){
              x<- "***"}else if(x <= 0.0001 && x > 0){
                x <- "****"}else if(x == 0){
                  x<- "****"}
  }))
  
  tukey_list$padj_significance2 <- unlist(lapply(tukey_list$padj, function(x) {
    
    if(x >= 0.1) {
      x <- ""} else if (x > 0.05 && x < 0.1){
        x <- format(x, digits=2, nsmall=2)} else if(x <= 0.05 && x > 0.01){
          x <- ""} else if (x <= 0.01 && x>0.001){
            x <- ""}else if ( x <= 0.001 && x>0.0001){
              x<- ""}else if(x <= 0.0001 && x > 0){
                x <- ""}else if(x == 0){
                  x<- ""}
  }))
  
  
  
  results <- list(anova=anovas, tukey = tukeys, tukey_list=tukey_list)
  return(results)
  
  
  
}
