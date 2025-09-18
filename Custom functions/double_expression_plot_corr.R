### double gene expression plots

double_expression_plot_corr <- function(gene_expression, gene_description,islet_id, genes_interest=c("INS", "IAPP"), folder=NULL,
                                   width=16, height=14, axis_size=16, legend_size=16, axis_ticks_size=14 ,
                                   x_title, color, ylim=c(0,15), xlim=c(0,15)){
  
  genes_interest <- data.frame(genesymbol = genes_interest)
  genes_list <- gene_description[gene_description$genesymbol %in% genes_interest$genesymbol,]
  genes_list <- genes_list[order(match(genes_list$genesymbol, genes_interest$genesymbol)),]
  
  if(nrow(genes_list) > length(genes_interest)){
    
    genes_list <- genes_list[!duplicated(genes_list[c("genesymbol")]),]
  }
  
  genes_interest_expression <- gene_expression[rownames(gene_expression) %in% rownames(genes_list),]
  genes_interest_expression <- genes_interest_expression[order(match(rownames(genes_interest_expression), rownames(genes_list))),]
  
  genes_interest_expression <- t(genes_interest_expression)
  genes_interest_expression <- as.data.frame(genes_interest_expression)
  colnames(genes_interest_expression) <- rownames(genes_list)
  genes_interest_expression$Sample <- rownames(genes_interest_expression)
  
  
  
  
  #longtable
  genes_interest_expression<- merge(genes_interest_expression, islet_id)
  
  
  ## plotting
  
  plot <- ggscatter (genes_interest_expression, 
                     x= colnames(genes_interest_expression)[2], 
                     y= colnames(genes_interest_expression)[3], 
                     fill = "Clinical_phenotype",
                     palette = color,
                     shape = "Clinical_phenotype",
                     size = 3,
                     #margin.plot = "boxplot",
                     add="reg.line", 
                     cor.coef = TRUE,
                     cor.coef.size = (axis_ticks_size/3 + 1),
                     ylab = genes_list$genesymbol[2],
                     xlab = genes_list$genesymbol[1],
                     # facet.by = "Marker",
                     # scales="free_x",
                     title = x_title,
                     font.title = c(legend_size, "bold", "darkblue"),
                     cor.method = "pearson",
                     add.params = list(color = "blue", size=1),
                     cor.coeff.args = list(method = "pearson", label.x.npc = "left", label.y.npc ="top"),
                     # panel.labs.font = list(face = "bold", color = "darkblue", size = 11, angle = NULL),
                     mean.point = FALSE, ellipse = FALSE, ellipse.type = "convex",
                     ggtheme = theme_bw() ) 
  
  plot+
    theme(axis.title.x = element_text(face = "bold.italic", size = axis_size), axis.title.y = element_text(face = "bold.italic", size = axis_size), 
          legend.title = element_text(face="bold", size = legend_size), legend.text = element_text(size = legend_size),
          axis.text = element_text(size = axis_ticks_size), plot.title = element_text(face = "bold", hjust = 0.5, size = legend_size) )+
    ylim(ylim)+ 
    #xlim(xlim)+
    scale_shape_manual(values=c(21:24)) 
  

  
}
