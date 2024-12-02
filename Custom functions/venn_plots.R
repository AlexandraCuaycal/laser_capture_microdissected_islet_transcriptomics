## custom venns
venn_plots <- function(data, DE_contrasts, folder_name, title,
                       adj.p.val = 0.001, logFCcutoff = 0.5,
                       width=15, height=13, dpi=300, units= "cm"){
  
  dir.create(paste(folder_name, "/Venn_plots", sep="" ))
  
  DEG_object <-data
  list_venn <- list()
  
  for (i in 1:length(DE_contrasts)){
    
    index_sig_DEGs <- which(DEG_object$resTable[[i]]$adj.P.Val < adj.p.val)
    index_logFC <- which(abs(DEG_object$resTable[[i]]$logFC) > logFCcutoff )
    index_DEGs <- intersect(index_sig_DEGs, index_logFC)
    sig_DEGs_table <- DEG_object$resTable[[i]][index_DEGs,]
    # ensembl_IDs <- ensembl_gene_id[rownames(ensembl_gene_id) %in% sig_DEGs_table$ProbesetID , ]
    # ensembl_IDs <-ensembl_IDs[ match(rownames(ensembl_IDs),sig_DEGs_table$ProbesetID ),]
    # sig_DEGs_table$ensembl_gene_ID <- ensembl_IDs$ensembl_gene_id
    list_venn[[i]] <- sig_DEGs_table[sort(sig_DEGs_table$adj.P.Val, index.return=TRUE)$ix,]
    list_venn[[i]] <- list_venn[[i]][order(list_venn[[i]]$logFC, decreasing = TRUE),]
    list_venn[[i]] <- list_venn[[i]]$ProbesetID
    
    
  }
  names(list_venn)<- DE_contrasts
  
  ggVennDiagram(x=list_venn, label_alpha = 0, category.names = change_to_groups(names(list_venn)),
                label = "both", set_size = 5, label_size = 5, label_percent_digit = 1) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    ggplot2::scale_fill_gradient(low="#8797E1",
                                 # mid="#DEDEE3", 
                                 high ="#E2889B"
                                 #,  midpoint = max(as.numeric(unlist(lapply(list_venn,length))))/2 
    )+
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title=element_text(size=13), 
          legend.text=element_text(size=11)) +
    labs(title = title )
  
  ggsave(filename = paste(folder_name,"/Venn_plots/Venn diagram_", title, ".tiff",sep = ""), width = width, height = height, dpi = dpi, units = units)
}