### double gene expression plots

double_expression_plot <- function(gene_expression, gene_description,islet_id, genes=c("INS", "IAPP"), folder=NULL,
                                   width=16, height=14, axis_size=16, legend_size=16, axis_ticks_size=14 ){
  
  genes_interest <- data.frame(probe_id = genes)
  genes_list <- gene_description[rownames(gene_description) %in% genes_interest$probe_id,]
  genes_list <- genes_list[order(match(genes_list$genesymbol, genes_interest$genesymbol)),]
  
  if(nrow(genes_list) > length(genes)){
    
    genes_list <- genes_list[!duplicated(genes_list[c("genesymbol")]),]
  }
  
  genes_interest_expression <- gene_expression[rownames(gene_expression) %in% rownames(genes_list),]
  genes_interest_expression <- genes_interest_expression[order(match(rownames(genes_interest_expression), rownames(genes_list))),]
  
  genes_interest_expression <- t(genes_interest_expression)
  genes_interest_expression <- as.data.frame(genes_interest_expression)
  colnames(genes_interest_expression) <- genes
  genes_interest_expression$Sample <- rownames(genes_interest_expression)
  
  
  
  
  #longtable
  genes_interest_expression<- merge(genes_interest_expression, islet_id)
  
  
  ## plotting
  
  plot <- ggplot(genes_interest_expression, aes(x=genes_interest_expression[,2], y=genes_interest_expression[,3], fill = Clinical_phenotype))
  plot + 
    geom_point(aes(shape=Clinical_phenotype), size=2)+
    #geom_points(aes( color =Clinical_phenotype, shape=Clinical_phenotype))+ 
    #geom_smooth(method = "lm", se=FALSE)+
    #scale_color_manual(values = color_pal)+
    theme_bw() +
    xlab(genes_list$genesymbol[1])+ylab(genes_list$genesymbol[2])+
    scale_shape_manual(values=c(21:24))+
    theme(axis.title.x = element_text(face = "bold.italic", size = axis_size), axis.title.y = element_text(face = "bold.italic", size = axis_size), 
          legend.title = element_text(face="bold", size = legend_size), legend.text = element_text(size = legend_size),
          axis.text = element_text(size = axis_ticks_size)) +
    scale_fill_manual(values=color_pal)
  
  
  
  
  
  
}