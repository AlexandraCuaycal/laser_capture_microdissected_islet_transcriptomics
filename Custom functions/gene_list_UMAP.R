## funtion to make the plots based on a gene list

gene_list_UMAP <- function( dataset, gene_list, gene_description_data=gene_description_data, colvec=color_pal, scale=2, controlscale=TRUE, new_folder=TRUE, folder_name = "folder", 
                            axistextsize = 10, legendtextsize = 10, dotsize=1.5,seed=175,
                            units = "cm", dpi = 300, width = 16, height = 12){
  
  genes_list_names <- gene_description_data[trimws(gene_description_data$genesymbol) %in% trimws(gene_list$genesymbol),]
  gene_probe_names <- rownames(genes_list_names)
  
  if(new_folder == TRUE){
    
    dir.create(folder_name)
    for(a in 1:length(gene_probe_names)){
      
      gene_interest <- as.vector(scale(as.numeric(exprs(dataset)[row.names(exprs(dataset)) == gene_probe_names[a],])))
      
      
      gene_umap <- umap(exprs(dataset), labels=gene_interest, scale=scale, controlscale=controlscale, colvec=colvec, dotsize=dotsize, seed=seed,
                        axistextsize = axistextsize, legendtextsize = legendtextsize, legendtitle = genes_list_names$genesymbol[a])
      
      ggsave(plot = gene_umap, filename = paste(folder_name, "/","UMAP_by_gene_", genes_list_names$genesymbol[a], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      
      print(paste("UMAP plot ", genes_list_names$genesymbol[a], " done", sep = ""))
      
    }
    
  } else{
    for(a in 1:length(gene_probe_names)){
      
      gene_interest <- as.vector(scale(as.numeric(exprs(dataset)[row.names(exprs(dataset)) == gene_probe_names[a],])))
      
      
      gene_umap <- umap(exprs(dataset), labels=gene_interest, scale=scale, controlscale=controlscale, colvec=colvec, dotsize=dotsize,
                        axistextsize = axistextsize, legendtextsize = legendtextsize, legendtitle = genes_list_names$genesymbol[a])
      
      ggsave(plot = gene_umap, filename = paste(folder_name, "/","UMAP_by_gene_", genes_list_names$genesymbol[a], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      
      print(paste("UMAP plot ", genes_list_names$genesymbol[a], " done", sep = ""))
      
    }
  }
  
  
}
