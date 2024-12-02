### lollipop plot

lollipop_plot <- function(result, title, label_size=11, showcat=10, cutstring=40, gene_number=40){
  
  df <- result@result
  
  df <- df %>% dplyr::filter(setSize >= gene_number)
  
  df$mycolor <- ifelse(df$NES<0, "cornflowerblue","indianred2")
  
  df <- df[order(df$p.adjust, decreasing = F),] #revert the order of the rows
  df$Label <-stringr::str_wrap(df$Description, width = cutstring) # keep max 50 character in a string (nchar("string") to count characters)
  df$Label <- factor(df$Label, levels = c(df$Label))
  
  df <- df %>% dplyr::group_by(mycolor) %>% slice_head(n = showcat)
  
  # Use limits of the x-axis range according to the range of p.values:
  max_x <- round(max(-log10(df$p.adjust)), digits = 0) + 0.5
  xlimits <- c(0, max_x)
  xbreaks <- seq_along(0:max_x)
  
  # Create the plot:
  p <- ggplot(df , aes(x = Label, y = -log10(p.adjust), fill = mycolor)) +
    geom_segment(aes(x = Label, xend = Label, y = 0, yend = -log10(p.adjust)),
                 color = df$mycolor, lwd = 1) +
    geom_point(pch = 21,  bg = df$mycolor, aes(size=df$setSize), color=df$mycolor) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_size(name = "Number\nof genes") + 
    scale_y_continuous(name=expression("-"*"log"[10]*"(adj p-value)"), 
                       limits=xlimits,
                       breaks=xbreaks) +
    scale_x_discrete(name="") +
    ggtitle(stringr::str_wrap(title, width = cutstring)) +
    theme_bw(base_size=10, base_family = "Helvetica") +
    theme(axis.text=element_text(size=label_size-1, colour = "black"),
          axis.title=element_text(size=label_size), legend.title=element_text(size=label_size-1), 
          legend.text=element_text(size=label_size-2), title = element_text(face = "bold", size = label_size)) +
    
    coord_flip()
  
  return(p)
  
  
}