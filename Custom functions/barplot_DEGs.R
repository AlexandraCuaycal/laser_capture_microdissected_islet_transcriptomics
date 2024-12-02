## to create bar plot of DE geens

barplot_DEGs <- function(DEG_tables, Contrast, Title, save=FALSE,
                         width=16, height=12, dpi=300){
  
  names(DEG_tables) <- Contrast
  
  count_table <- NULL
  
  for(i in 1:length(DEG_tables) ) {
    
    DEG_tables[[i]] <- DEG_tables[[i]] %>% mutate(diffexpressed = case_when(
      logFC > 0  ~ 'UP',
      logFC < 0  ~ 'DOWN'))
    
    counts_g <- data.frame(Contrast=rep(Contrast[i], times=nrow(DEG_tables[[i]])), Diffexpressed = DEG_tables[[i]]$diffexpressed)
    count_table <- dplyr::bind_rows( count_table, counts_g )
    
  }
  
  count_table$Contrast <- factor(count_table$Contrast, levels = Contrast)
  count_table$Diffexpressed <- factor(count_table$Diffexpressed, levels=c("UP", "DOWN"))
  
  p <- ggplot(count_table, aes(x = Contrast, fill=Diffexpressed))
  
  p<- p + geom_bar() +theme_bw() +
    
    # geom_text(stat = "count", aes(label = after_stat(count)), position = position_stack(vjust = 0.5),
    #                 fontface="bold", size = 5)+
    
    geom_text_repel(stat = "count", aes(label = after_stat(count)), position = position_stack(vjust = 0.4),
                    fontface="bold", size = 5, direction = "y", force = 0.15, force_pull = 1.95,
                    max.iter = 100, max.time = 0.1)+
    
    ggtitle(Title)+
    theme(axis.text.x=element_text(size=14, face = "bold", angle = 25, hjust=1), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title=element_text(size=14,face="bold"), legend.text = element_text(size = 14),
          legend.title =element_text(size=14,face="bold"),
          plot.title = element_text(size = 14, face = "bold", hjust=0.5)) +
    scale_fill_brewer(palette = "Pastel1")
  
  plot(p)
  
  if(save){
    ggsave2(filename =paste(Title, ".tiff", sep = ""), width = width, height = height, dpi = dpi, units = "cm")
  }
  
  
}


