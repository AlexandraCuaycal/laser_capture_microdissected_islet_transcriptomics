
## these functions will make the plots with categorical variables as provided

## for PCA

pca_plots <- function(dataset, cat_var, shape_var, col_pal = color_pal, frame=FALSE, new_folder=TRUE, folder_name = "folder", units = "cm", 
                      dpi = 300, width = 16, height = 12 , scale=FALSE, loadings=FALSE,
                      loadings.label=FALSE, scale2=0,  frame.type='norm'){
  
  index_cat_var <- match(cat_var, varLabels(dataset))
  pdata <- pData(phenoData(dataset))[,index_cat_var]
  
  if(new_folder == TRUE){
    dir.create(folder_name)
    pca <- prcomp(t(exprs(dataset)), scale. = scale, center = TRUE)
    for(i in 1:length(cat_var)){
      
      autoplot(pca, data=pdata,  shape=shape_var,colour=cat_var[i],
               frame = frame,  frame.type =frame.type, frame.alpha=0.05,
               loadings=loadings, loadings.label = loadings.label,
               scale=scale2) +theme_bw() + 
        scale_color_manual(values = col_pal) + 
        scale_fill_manual(values = col_pal) 
      ggsave(filename = paste(folder_name, "/","PCA_by_", cat_var[i], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      print(paste("PCA ", cat_var[i], " done", sep = ""))
    }
  }else{
    pca <- prcomp(t(exprs(dataset)), scale. =  scale, center = TRUE)
    for(i in 1:length(cat_var)){
      
      autoplot(pca, data=pdata,  shape=shape_var,colour=cat_var[i],
               frame = frame,  frame.type = frame.type, frame.alpha=0.05,
               loadings=loadings, loadings.label = loadings.label,
               scale=scale2) +theme_bw() + 
        scale_color_manual(values = col_pal) + 
        scale_fill_manual(values = col_pal)  
      ggsave(filename = paste("PCA_by_", cat_var[i], ".tiff", sep = ""), width =width, height = height, units = units, dpi = dpi )
      print(paste("PCA ", cat_var[i], " done", sep = ""))
    }
  }
}