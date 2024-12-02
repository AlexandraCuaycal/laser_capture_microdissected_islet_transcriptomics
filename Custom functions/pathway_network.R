## pathway network plot

pathway_network <- function(result, title, label_size=11, showcat=10, gene_number=40, repel = F,
                            enrichment=F){
  
  if(enrichment == T){
    
    df <- result
    df <- df %>% dplyr::filter(Count >= gene_number)
    df$mycolor <- rep("lightblue", times = nrow(df))
    
  } else {
    df <- result@result
    df <- df %>% dplyr::filter(setSize >= gene_number)
    
    df$mycolor <- ifelse(df$NES<0, "cornflowerblue","indianred2")
  }
  
  
  
  
  df <- df[order(df$p.adjust, decreasing = F),] #revert the order of the rows
  
  df <- df %>% dplyr::group_by(mycolor) %>% slice_head(n = showcat)
  
  settings <- aPEAR.theme
  
  settings$fontSize <- label_size
  settings$repelLabels <- repel
  
  enrichmentNetwork(df, theme = settings) 
  
}
