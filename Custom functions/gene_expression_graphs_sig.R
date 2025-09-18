###### gene graphs

gene_expression_graphs_sig <- function(gene_expression, gene_description, genes_interest, folder=NULL, width=16, height=14,
                                   x_title="INS+ CD3- islets", color, islet_id, legend_size=16,
                                   axis_ticks_size=14, ylim=c(0, 16)) {
  
  if(length(genes_interest) == 1){
    
    genes_list <- gene_description_data[trimws(gene_description_data$genesymbol) == genes_interest,]
    
  }else{
    
    genes_list <- gene_description[rownames(gene_description) %in% genes_interest$probe_id,]
  }
  genes_list <- subset(genes_list, locus.type == "Coding", select= - start)
  genes_list$probe_id <- rownames(genes_list)
  
  
  # table with average and se
  genes_interest_expression <- table_ready(gene_expression = gene_expression, genes_of_interest =  genes_list)
  
  
  # table with expression from all islets
  
  if (nrow(genes_list) == 1){
    genes_interest_expression2 <- gene_expression[rownames(gene_expression) %in% rownames(genes_list),]
    
    genes_interest_expression2 <- as.data.frame(t(genes_interest_expression2))
    rownames(genes_interest_expression2) <- rownames(genes_list)
    genes_interest_expression2 <- merge(genes_interest_expression2, genes_list)
    #long table
    genes_interest_expression2 <- reshape2::melt(genes_interest_expression2, value.name = "Expression", variable.name = "Sample")
    
  }else {
    genes_interest_expression2 <- gene_expression[rownames(gene_expression) %in% rownames(genes_list),]
    genes_interest_expression2 <- genes_interest_expression2[order(match(rownames(genes_interest_expression2), rownames(genes_list))),]
    
    
    rownames(genes_interest_expression2) <- rownames(genes_list)
    genes_interest_expression2 <- as.data.frame(t(genes_interest_expression2))
    genes_interest_expression2$Sample <- rownames(genes_interest_expression2)
    
    # long expression table
    genes_interest_expression2 <- reshape2::melt(genes_interest_expression2, value.name="Expression", variable_name = "probe_id", 
                                                 id.vars = "Sample")
    colnames(genes_interest_expression2)[2]<- "probe_id"
    

    genes_interest_expression2 <- merge(genes_interest_expression2, genes_list)
  }
  

  #merge with islet metdata
  genes_interest_expression2<- merge(genes_interest_expression2, islet_id)
  
  ## stat test
  
  stat.test <- tukey_tests(gene_expression, genes_list)$tukey_list
  
  stat.test <- stat.test %>% filter(p.adj.signif != "ns")
  
  
  dat <- genes_interest_expression[trimws(genes_interest_expression$genesymbol) == genes_interest,]
  
  
  plot <- ggplot(data=dat)
  
  plot <- plot +
    geom_errorbar(data=dat, color="black", aes(x=Clinical_phenotype, ymin=Expression-Std_error, 
                                               ymax=Expression+Std_error, group=Clinical_phenotype), width= 0.06, position=position_dodge(0.05)) +
    geom_line(data=dat, color="gray29", aes(x=Clinical_phenotype, y=Expression, group=probe_id))+
    geom_point(data = genes_interest_expression2, aes(x=Clinical_phenotype, y = Expression),
               position = position_identity(), color='darkgrey', size=1) +
    
    geom_point(data=dat,aes(x=Clinical_phenotype, y=Expression, fill=Clinical_phenotype, group=probe_id, shape=Clinical_phenotype),size=4) + 
    
    
    ggpubr::stat_pvalue_manual(data=stat.test, label = "p.adj.signif", y.position = stat.test$y.position+0.05,
                               # tip.length=0.01
                               tip.length =   0.02, label.size = 4.5, step.increase = 0.05 )+
    
    theme_bw()+  theme( axis.text = element_text(size = axis_ticks_size, face = "bold"), axis.title = element_text(size=legend_size, face = "bold"), 
                        legend.text = element_text(size = legend_size), panel.grid.minor = element_blank(), 
                        legend.title = element_text(size=legend_size, face="bold"),
                        #axis.text.x = element_text( angle = 60, vjust = 0.5),
                        axis.title.x = element_text(color="darkblue")
                        
    ) + 
    ggtitle(genes_interest) + theme(plot.title = element_text(hjust = 0.5, size = legend_size, face = "bold.italic"))+
    
    xlab(paste("\n", x_title, sep = ""))+ scale_color_manual(values= color_pal
                                                             #, labels= c("ND (n=10)", "sAAb (n=3)", "mAAb (n=3)", "T1D (n=6)")
    )+
    
    scale_x_discrete(labels= stringr::str_wrap(c("ND", "sAAb", "mAAb", "T1D"), width = 5), 
                     #guide = ggplot2::guide_axis(n.dodge = 2)
                     
    )+
    
    #coord_cartesian(ylim=ylim)+

    scale_shape_manual(values=c(21:24))+
    scale_fill_manual(values=color_pal)
  
  
  
}
