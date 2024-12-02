### bar plot of -log10 p value

barplot_result <- function(result, title, color, showcat=15){
  
  df<- result@result
  
  
  df %>%
    dplyr::arrange(p.adjust) %>%
    slice_head(n = showcat) %>%
    ggplot(aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust))) +
    
    geom_bar(stat = "identity", fill =color) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    labs(y = "Description", x = expression("-"*"log"[10]*"(adj p-value)")) +
    ggtitle(title) +
    theme_classic()
  
}
