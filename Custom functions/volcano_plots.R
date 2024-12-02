#volcano plots

volcano_plots <- function(data, DE_contrasts, folder_name, title=title, labsize=3,
                          x="logFC", y="adj.P.Val", pCutoff = 1e-03, FCcutoff = 1.5,
                          width=20, height=22, dpi=300, units= "cm"){
  
  dir.create(paste(folder_name, "/Volcano_plots", sep="" ))
  
  
  for (i in 1:length(DE_contrasts)){
    data_for_plot <- data[["resTable"]][[DE_contrasts[i]]]
    ensembl_IDs <- ensembl_gene_id[rownames(ensembl_gene_id) %in% data_for_plot$ProbesetID , ]
    ensembl_IDs <-ensembl_IDs[ match(rownames(ensembl_IDs),data_for_plot$ProbesetID ),]
    data_for_plot$ensembl_gene_ID <- ensembl_IDs$ensembl_gene_id
    data_for_plot <- data_for_plot[grepl("ENSG", data_for_plot$ensembl_gene_ID),]
    
    EnhancedVolcano(data_for_plot, labSize = labsize,
                    lab = data_for_plot$GeneName, x=x, y=y, pCutoff = pCutoff, FCcutoff = FCcutoff,
                    title=title, subtitle = DE_contrasts[i])
    ggsave(filename = paste(folder_name, "/Volcano_plots/", DE_contrasts[i], ".tiff", sep="" ),
           width = width, height = height, dpi = dpi, units = units)
    
    print(paste(DE_contrasts, " done"))
    
  }
  
}
