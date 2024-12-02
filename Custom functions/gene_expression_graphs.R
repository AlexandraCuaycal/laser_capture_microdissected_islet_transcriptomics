###### gene graphs

gene_expression_graphs <- function(gene_expression, gene_description, genes_interest, folder=NULL, width=16, height=14,
                                   x_title="INS+ CD3- islets") {
  
  if(ncol(genes_interest) == 1){
    
    genes_list <- gene_description_data[trimws(gene_description_data$genesymbol) %in% genes_interest$genesymbol,]
    genes_list <- subset(genes_list, locus.type == "Coding")
  }else{
    
    genes_list <- gene_description[rownames(gene_description) %in% genes_interest$probe_id,]
  }
  
  
  genes_interest_expression <- table_ready(gene_expression = gene_expression, genes_of_interest =  genes_list)
  
  tukey_list <- tukey_tests(gene_expression, genes_list)$tukey_list
  
  if(!is.null(folder)){
    
    dir.create(folder)
    
    for (m in 1:length(genes_interest$genesymbol)){
      
      plot <- ggplot(genes_interest_expression[trimws(genes_interest_expression$genesymbol) == genes_interest$genesymbol[m],], 
                     aes(x=Clinical_phenotype, y=Expression, fill=Clinical_phenotype, group=probe_id, label=genesymbol))
      
      ## for sig stars
      significance <- tukey_list [ trimws(tukey_list$genesymbol) == genes_interest$genesymbol[m],]
      significance <- significance[unlist(lapply(significance$Comparison, function(x) {
        grepl("ND", x, fixed = TRUE)
      })),]
      
      ## for p value of almost sig
      significance2 <- tukey_list [ trimws(tukey_list$genesymbol) == genes_interest$genesymbol[m],]
      significance2 <- significance2[unlist(lapply(significance2$Comparison, function(x) {
        grepl("ND", x, fixed = TRUE)
      })),]
      
      
      significance$Comparison <- factor(significance$Comparison, levels = c("sAAb-ND", "ND-mAAb" , "T1D-ND"))
      significance <- with(significance, significance[order(Comparison),])
      
      significance2$Comparison <- factor(significance2$Comparison, levels = c("sAAb-ND", "ND-mAAb" , "T1D-ND"))
      significance2 <- with(significance2, significance2[order(Comparison),])
      
      
      
      plot +
        theme_bw()+  theme( axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size=12, face = "bold"), 
                            legend.text = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(), 
                            legend.title = element_text(size=12, face="bold"),
                            #axis.text.x = element_text(angle=60, hjust=1)
        ) + 
        
        xlab(paste("\n", x_title, sep = ""))+ scale_color_manual(values= color_pal
                                                                 #, labels= c("ND (n=10)", "sAAb (n=3)", "mAAb (n=3)", "T1D (n=6)")
        )+ 
        
        scale_x_discrete(labels= stringr::str_wrap(c("ND", "sAAb", "mAAb", "T1D"), width = 5), 
                         #guide = ggplot2::guide_axis(n.dodge = 2)
                         
        )+
        
        #geom_errorbar(color="darkgrey", aes( ymin=Expression-Std_error, ymax=Expression+Std_error, group=Clinical_phenotype), width= 0.06, position=position_dodge(0.05)) +
        geom_line(color="darkgrey")+
        geom_point(size=4, aes(shape=Clinical_phenotype)) + ggtitle(genes_interest$genesymbol[m]) + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold.italic"))+
        
        scale_shape_manual(values=c(21:24))+
        scale_fill_manual(values=color_pal)+
        # gene symbol
        #geom_text( aes(x=Clinical_phenotype, y=Expression+Expression/1000, label=ifelse(Expression >= 0, as.character(genesymbol), ''), hjust=0.5, vjust=-0.5), size = 4, fontface = "italic",show.legend = FALSE) +
        
        ## adding significance to plots
        geom_text( aes(x=Clinical_phenotype, y=Expression +Expression/1000, label= ifelse(Expression != 0,c(rep("", times = length(unique(significance$probe_id))), significance$padj_significance), ""), hjust=0.5, vjust=-0.05), colour= "#282828", fontface="bold",  size = 7, show.legend = FALSE) +
        
        ## adding pvalues for almost significance
        geom_text( aes(x=Clinical_phenotype, y=Expression +Expression/1000, label= ifelse(Expression != 0,c(rep("", times = length(unique(significance2$probe_id))), significance2$padj_significance2), ""), hjust=0.5, vjust=-1), colour= "#282828", fontface="bold",  size = 4, show.legend = FALSE)
      
      
      #geom_text( aes(x=Clinical_phenotype, y=Expression, label=ifelse(Expression >= 0, as.character(probe_id), ''), hjust=0, vjust=-1), size = 1.5, fontface = "italic",show.legend = FALSE)
      print(genes_interest$genesymbol[m])
      #ggsave(filename = paste("Mollie - gene expression graphs/Not insulitic/",genes_interest$genesymbol[m], "_gene expresion_graph.tiff", sep = ""), height = 14, width = 16, units = "cm", dpi=300)
      ggsave(filename = paste(genes_interest$genesymbol[m], "_gene expresion_graph.tiff", sep = ""), height = height, width = width, units = "cm", dpi=300)
      
      
    }
    
    
  }else{
    for (m in 1:length(genes_interest$genesymbol)){
      
      plot <- ggplot(genes_interest_expression[trimws(genes_interest_expression$genesymbol) == genes_interest$genesymbol[m],], 
                     aes(x=Clinical_phenotype, y=Expression, fill=Clinical_phenotype, group=probe_id, label=genesymbol))
      
      ## for sig stars
      significance <- tukey_list [ trimws(tukey_list$genesymbol) == genes_interest$genesymbol[m],]
      significance <- significance[unlist(lapply(significance$Comparison, function(x) {
        grepl("ND", x, fixed = TRUE)
      })),]
      
      ## for p value of almost sig
      significance2 <- tukey_list [ trimws(tukey_list$genesymbol) == genes_interest$genesymbol[m],]
      significance2 <- significance2[unlist(lapply(significance2$Comparison, function(x) {
        grepl("ND", x, fixed = TRUE)
      })),]
      
      
      significance$Comparison <- factor(significance$Comparison, levels = c("sAAb-ND", "ND-mAAb" , "T1D-ND"))
      significance <- with(significance, significance[order(Comparison),])
      
      significance2$Comparison <- factor(significance2$Comparison, levels = c("sAAb-ND", "ND-mAAb" , "T1D-ND"))
      significance2 <- with(significance2, significance2[order(Comparison),])
      
      
      
      plot<- plot +
        theme_bw()+  theme( axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size=12, face = "bold"), 
                            legend.text = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(), 
                            legend.title = element_text(size=12, face="bold"),
                            #axis.text.x = element_text(angle=60, hjust=1)
        ) + 
        
        xlab(paste("\n", x_title, sep = ""))+ scale_color_manual(values= color_pal
                                                                 #, labels= c("ND (n=10)", "sAAb (n=3)", "mAAb (n=3)", "T1D (n=6)")
        )+ 
        
        scale_x_discrete(labels= stringr::str_wrap(c("ND", "sAAb", "mAAb", "T1D"), width = 5), 
                         #guide = ggplot2::guide_axis(n.dodge = 2)
                         
        )+
        
        #geom_errorbar(color="darkgrey", aes( ymin=Expression-Std_error, ymax=Expression+Std_error, group=Clinical_phenotype), width= 0.06, position=position_dodge(0.05)) +
        geom_line(color="darkgrey")+
        geom_point(size=4, aes(shape=Clinical_phenotype)) + ggtitle(genes_interest$genesymbol[m]) + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold.italic"))+
        
        # gene symbol
        #geom_text( aes(x=Clinical_phenotype, y=Expression+Expression/1000, label=ifelse(Expression >= 0, as.character(genesymbol), ''), hjust=0.5, vjust=-0.5), size = 4, fontface = "italic",show.legend = FALSE) +
        
        scale_shape_manual(values=c(21:24))+
        scale_fill_manual(values=color_pal)+
        ## adding significance to plots
        geom_text( aes(x=Clinical_phenotype, y=Expression +Expression/1000, label= ifelse(Expression != 0,c(rep("", times = length(unique(significance$probe_id))), significance$padj_significance), ""), hjust=0.5, vjust=-0.05), colour= "#282828", fontface="bold",  size = 7, show.legend = FALSE) +
        
        ## adding pvalues for almost significance
        geom_text( aes(x=Clinical_phenotype, y=Expression +Expression/1000, label= ifelse(Expression != 0,c(rep("", times = length(unique(significance2$probe_id))), significance2$padj_significance2), ""), hjust=0.5, vjust=-1), colour= "#282828", fontface="bold",  size = 4, show.legend = FALSE)
      
      
      #geom_text( aes(x=Clinical_phenotype, y=Expression, label=ifelse(Expression >= 0, as.character(probe_id), ''), hjust=0, vjust=-1), size = 1.5, fontface = "italic",show.legend = FALSE)
      print(genes_interest$genesymbol[m])
      
      plot(plot)
      #ggsave(filename = paste("Mollie - gene expression graphs/Not insulitic/",genes_interest$genesymbol[m], "_gene expresion_graph.tiff", sep = ""), height = 14, width = 16, units = "cm", dpi=300)
      #ggsave(filename = paste(genes_interest$genesymbol[m], "_gene expresion_graph.tiff", sep = ""), height = height, width = width, units = "cm", dpi=300)
      
      
    }
  }
  
  
  
}