### PCA or UMAP plots in M3C package
#23, 27,37
PCA_UMAP_plots <- function(dataset, cat_var, colvec=color_pal, seed=37,scale=3, controlscale=TRUE, new_folder=TRUE, folder_name = "folder", 
                           axistextsize = 10, legendtextsize = 10, dotsize = 1.5,
                           units = "cm", dpi = 300, width = 16, height = 12,
                           scaler=FALSE){
  
  index_cat_var <- match(cat_var, varLabels(dataset))
  pdata <- pData(phenoData(dataset))[,index_cat_var]
  pdata <- as.data.frame(pdata)
  
  if(new_folder == TRUE){
    dir.create(folder_name)
    for(i in 1:length(cat_var)){
      
      pca_p <- M3C::pca(exprs(dataset), labels=pdata[,i],scale=scale, controlscale=controlscale, colvec=colvec, dotsize=dotsize,
                   axistextsize = axistextsize, legendtextsize = legendtextsize, legendtitle = colnames(pdata)[i],
                   scaler = scaler)

      umap_p <- M3C::umap(exprs(dataset), seed=seed,labels=pdata[,i],scale=scale, controlscale=controlscale, colvec=colvec, dotsize=dotsize,
                     axistextsize = axistextsize, legendtextsize = legendtextsize, legendtitle = colnames(pdata)[i])
      
      plot(pca_p)
      ggsave(plot = pca_p, filename = paste(folder_name, "/","PCA_by_", cat_var[i], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      
     # pdf(file = paste(folder_name, "/","UMAP_by_", cat_var[i], ".pdf", sep = ""), width =width, height = height)
      plot(umap_p)
      ggsave(plot = umap_p, filename = paste(folder_name, "/","UMAP_by_", cat_var[i], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      #dev.off()
      print(paste("Plots ", cat_var[i], " done", sep = ""))
      
    }
    
  }else{
    for(i in 1:length(cat_var)){
      
      pca_p <- M3C::pca(exprs(dataset), labels=pdata[,i],scale=scale, controlscale=controlscale, colvec=colvec, dotsize=dotsize,
                   axistextsize = axistextsize, legendtextsize = legendtextsize, legendtitle = colnames(pdata)[i],
                   scaler = scaler)

      umap_p <- M3C::umap(exprs(dataset), seed=seed, labels=pdata[,i],scale=scale, controlscale=controlscale, colvec=colvec, dotsize=dotsize,
                     axistextsize = axistextsize, legendtextsize = legendtextsize, legendtitle = colnames(pdata)[i])
      
      plot(pca_p)
      ggsave(plot = pca_p, filename = paste(folder_name, "/","PCA_by_", cat_var[i], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      
      plot(umap_p)
      ggsave(plot = umap_p, filename = paste(folder_name, "/","UMAP_by_", cat_var[i], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      
      print(paste("Plots ", cat_var[i], " done", sep = ""))
      
    }
  }
  
  
}