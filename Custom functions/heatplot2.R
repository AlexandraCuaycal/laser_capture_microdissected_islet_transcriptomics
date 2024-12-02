

heatplot.enrichResult <- function(x, showCategory=30, foldChange=NULL,
                                  label_format = 30) {
  
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  n <- update_n(x, showCategory)
  geneSets <- extract_geneSets(x, n)
  
  foldChange <- fc_readable(x, foldChange)
  d <- list2df(geneSets)
  
  if (!is.null(foldChange)) {
    d$foldChange <- foldChange[as.character(d[,2])]
    ## palette <- fc_palette(d$foldChange)
    p <- ggplot(d, aes_( ~categoryID, ~Gene)) +
      geom_tile(aes_(fill = ~foldChange), color = "white") +
      xlab("")+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 15))+
      scale_fill_continuous(low="blue", high="red", name = "fold change")
    ## scale_fill_gradientn(name = "fold change", colors = palette)
    
  } else {
    p <- ggplot(d, aes_(~categoryID, ~Gene)) + geom_tile(color = 'white')
  }
  p + xlab(NULL) + ylab(NULL) + theme_minimal() +
    scale_y_discrete(labels = label_func) +
    xlab("")+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15))+
    theme(panel.grid.major = element_blank(),
          axis.text.x=element_text(angle = 60, hjust = 1))
}