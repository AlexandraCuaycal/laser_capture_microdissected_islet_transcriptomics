lollipop_plot_enrich <- function(result, title, color, showcat=15, cutstring = 40,
                                 label_size = 11){
  
  df<- result@result
  
  df$Label <-stringr::str_trunc(df$Description, cutstring, "right") # keep max 50 character in a string (nchar("string") to count characters)
  df$Label <- factor(df$Label, levels = c(df$Label))
  
  # Use limits of the x-axis range according to the range of p.values:
  max_x <- round(max(-log10(df$p.adjust)), digits = 0) + 0.5
  xlimits <- c(0, max_x)
  xbreaks <- seq_along(0:max_x)
  
  
  # Create the plot:
  p <- ggplot(df %>%
                dplyr::arrange(p.adjust) %>%
                slice_head(n = showcat) , aes(x = Label, y = -log10(p.adjust), colour = color)) +
    geom_segment(aes(x = Label, xend = Label, y = 0, yend = -log10(p.adjust)),
                 color = color, lwd = 1) +
    geom_point(pch = 21,  bg = color, aes(size=Count), color=color) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_size(name = "Number\nof genes") + 
    scale_y_continuous(name=expression("-"*"log"[10]*"(adj p-value)"), 
                       limits=xlimits,
                       breaks=xbreaks) +
    scale_x_discrete(name="") +
    theme_bw(base_size=10, base_family = "Helvetica") +
    theme(axis.text=element_text(size=label_size, colour = "black"),
          axis.title=element_text(size=label_size), legend.title=element_text(size=label_size-1), 
          legend.text=element_text(size=label_size-2), title = element_text(face = "bold")) +
    ggtitle(title) +
    coord_flip()
  
  return(p)
  
}