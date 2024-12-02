barplot_nes <- function(result, title, showcat){
  
  sorted_result<- result@result[order(result@result$NES, decreasing = F),]
  
  sorted_result$color<-ifelse(sorted_result$NES<0, "indianred2", "cornflowerblue")
  
  sorted_result %>%
    dplyr::group_by(color) %>%
    dplyr::arrange(desc(abs(NES))) %>%
    slice_head(n = showcat) %>%
    ggplot(aes(x = NES, y = reorder(Description, NES), fill = color)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0) +
    labs(y = "Description") +
    ggtitle(title)+
    theme_classic() +
    theme(legend.position = "none")
}