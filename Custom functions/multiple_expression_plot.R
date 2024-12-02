multiple_expression_plot <- function(gene_expression, gene_description,islet_id, genes_interest, folder=NULL,
                                     width=16, height=14, label_size=5, legend_size=14, title,
                                     title_size=16, angle=0){
  
  genes_list <- gene_description[rownames(gene_description) %in% genes_interest$probe_id,]
  genes_list <- genes_list[order(match(genes_list$genesymbol, genes_interest$genesymbol)),]
  
  if(nrow(genes_list) > length(genes_interest$genesymbol)){
    
    genes_list <- genes_list[!duplicated(genes_list[c("genesymbol")]),]
  }
  
  genes_interest_expression <- gene_expression[rownames(gene_expression) %in% rownames(genes_list),]
  genes_interest_expression <- genes_interest_expression[order(match(rownames(genes_interest_expression), rownames(genes_list))),]
  
  genes_interest_expression <- as.data.frame(genes_interest_expression)
  genes_interest_expression$Probe_id <- rownames(genes_interest_expression)
  genes_interest_expression$genesymbol <- genes_list$genesymbol
  
  
  #longtable
  genes_interest_expression_long <-  reshape2::melt(genes_interest_expression, value.name =  "Expression", variable.name="Sample")
  
  genes_interest_expression_long<- merge(genes_interest_expression_long, islet_id)
  
  genes_interest_expression_long <- genes_interest_expression_long %>% group_by(Clinical_phenotype, genesymbol) %>%
    mutate(Avg_expression = mean(Expression))
  
  genes_interest_expression_long$Size <- scale_values(genes_interest_expression_long$Avg_expression)
  
  genes_interest_expression_long$Clinical_phenotype<-factor(genes_interest_expression_long$Clinical_phenotype, 
                                                            levels=c("ND", "sAAb", "mAAb", "T1D"))
  genes_interest_expression_long$genesymbol <- factor(genes_interest_expression_long$genesymbol, 
                                                      levels = genes_list$genesymbol)
  
  ## for polar plot
  
  data_radar <- genes_interest_expression_long %>% group_by(Clinical_phenotype, Sample) %>% dplyr::select(Avg_expression, genesymbol) %>%
    pivot_wider(values_from =Avg_expression, names_from = c(genesymbol))
  
  data_radar <- subset(data_radar[!duplicated(data_radar[,"Clinical_phenotype"]),], select=-Sample)
  
  data_radar_scaled <- data_radar %>% as_tibble(rownames = "Clinical_phenotype") %>% 
    mutate(across(-Clinical_phenotype, scales::rescale 
                  #, from = c(0,max(apply(data_radar[,-1], 2, max)))
                  
    ))
  
  
  
  r<-ggradar_edit(data_radar_scaled, group.colours = color_pal, group.line.width = 1, group.point.size = 4,
                  legend.title = "Clinical_phenotype", legend.position = "bottom", axis.label.font='italic',
                  values.radar = rep("", times=3), axis.label.size = label_size , angle=angle)
  
  # r<-ggradar_edit(data_radar, group.colours = color_pal, group.line.width = 1, group.point.size = 4,
  #                 legend.title = "Clinical_phenotype", legend.position = "bottom", axis.label.font='italic')  
  # 
  
  r<- r + ggtitle(title ) +
    theme(legend.title = element_text(size=legend_size), legend.text = element_text(size=legend_size),
                plot.title = element_text(hjust = 0.5, size = title_size, face = "bold"))
    
    
  
  ## plotting
  
  # p <- ggplot(genes_interest_expression_long , aes(x=genesymbol, y=Clinical_phenotype))
  # 
  # p<- p+
  #   geom_point(aes(size=Size, color=Avg_expression)) +theme_bw()+theme(axis.title = element_blank(), axis.text.x = element_text(angle = 60, vjust = 0.5),
  #                                                                    axis.text = element_text(size = 12),
  #                                                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                                                                    panel.background = element_blank())+
  #   scale_color_gradient(low = "lightgrey", high = "blue", name = "Avg expression",
  #                          limits=c(0,max(genes_interest_expression_long$Avg_expression)+0.5)
  #                          )+
  #   scale_size_continuous(guide = NULL)
  # 
  # plot(p)
  
  plot(r)
  
  
}
