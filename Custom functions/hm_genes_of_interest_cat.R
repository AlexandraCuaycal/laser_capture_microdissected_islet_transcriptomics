hm_genes_of_interest_cat <- function(genes_of_interest, gene_description_data, gene_expression_obj, 
                                 group_of_interest = "sAAb", Insulitis_cat = "Not_insulitic", histology=TRUE, colors=color_pal,
                                 width=20, height=10, process="Complement", by_donor=TRUE, R=FALSE, res=600,
                                 C=FALSE,
                                 min=-2, max=2, rowname=10, rowtitle=11){
  
  
  ## getting only the INS pos islets
  
  gene_expression_obj <- gene_expression_obj[,gene_expression_obj$INS_stain == "Pos"]
  
  
  
  gene_expression <- exprs(gene_expression_obj)
  
  # metdata
  islet_id <- pData(phenoData(gene_expression_obj))
  
  islet_id$Islet_N <- rownames(islet_id)
  
  islet_id<- islet_id %>% group_by(Clinical_phenotype) 
  
  islet_id <- islet_id[order(islet_id$Clinical_phenotype),]
  
  rownames(islet_id) <- islet_id$Islet_N
  
  gene_expression<- gene_expression[, match(islet_id$Islet_N, colnames(gene_expression))]
  
  # for donor type of interest
  islet_id_interest<- islet_id %>% filter(Clinical_phenotype %in% c("ND",group_of_interest))
  
  if (histology == TRUE){
    group_histo  <- islet_id_interest %>% filter(Clinical_phenotype %in% group_of_interest, 
                                                 Insulitis == Insulitis_cat) %>%
      dplyr::select(Islet_N)
  }else{
    group_histo  <- islet_id_interest %>% filter(Clinical_phenotype %in% group_of_interest) %>%
      dplyr::select(Islet_N)
  }
  
  
  
  group_ND <- islet_id_interest %>% filter(Clinical_phenotype == "ND", 
                                           Insulitis == "Not_insulitic") %>%
    dplyr::select(Islet_N)
  
  islet_id_interest <- islet_id_interest %>% filter(Islet_N %in% c(group_ND$Islet_N, group_histo$Islet_N))
  islet_id_interest <- islet_id_interest %>% group_by(Clinical_phenotype) %>% arrange(Histological_phenotype)
  
  
  gene_expression <- gene_expression[, colnames(gene_expression) %in% islet_id_interest$Islet_N]
  gene_expression<- gene_expression[, match(islet_id_interest$Islet_N, colnames(gene_expression))]
  
  
  expression_list <- gene_expression[rownames(gene_expression) %in% genes_of_interest$probe_id,]
  
  expression_list<- expression_list[match(genes_of_interest$probe_id, rownames(expression_list)),]
  
  expression_list_cat <- data.frame(expression_list, genesymbol=genes_of_interest$genesymbol,
                                    Category = genes_of_interest$Category)
  
  ## calculating donor averages
  
  expression_list_donor <- NULL
  donors <- unique(islet_id_interest$Donor_ID)
  for(i in 1:length(donors)){
    
    islets_donor <- islet_id_interest %>% filter(Donor_ID == donors[i]) %>%
      dplyr::select(Islet_N)
    
    expression_islets <- expression_list[, colnames(expression_list) %in% islets_donor$Islet_N]
    
    if(! is.vector(expression_islets) == T){
      
      expression_islets <- apply(expression_islets, 1, mean)
    }
    
    
    
    expression_list_donor <- cbind(expression_list_donor, expression_islets)
    
    
  }
  colnames(expression_list_donor)<- donors
  
  
  
  expression_list_donor_cat <- data.frame(expression_list_donor, genesymbol=genes_of_interest$genesymbol,
                                      Category = genes_of_interest$Category)
  

  
  
  
  
  ## heatmap
  library(RColorBrewer)
  require(ComplexHeatmap)
  require(circlize)
  require(digest)
  require(cluster)
  
  
  ## we need factor variables: categories, colnames same
  
    # to calculate z scores
  
  if (by_donor == TRUE){
    x <- expression_list_donor_cat[,(1:ncol(expression_list_donor))]
    
    
    islet_id_interest <- islet_id_interest[!duplicated(islet_id_interest$Donor_ID),]
    
  }else{
    x <- expression_list_cat[,(1:ncol(expression_list))]
  }
  
  
  
  nsamples <- ncol(x)
  
  M <- rowMeans(x, na.rm = TRUE)
  DF <- nsamples - 1L
  IsNA <- is.na(x)
  if (any(IsNA)) {
    mode(IsNA) <- "integer"
    DF <- DF - rowSums(IsNA)
    DF[DF == 0L] <- 1L
  }
  x <- x - M
  V <- rowSums(x^2L, na.rm = TRUE)/DF
  x <- x/sqrt(V + 0.01)
  
  
  #save(x, file = "z_scores_data.Rdata")
  
  
  ####
  
  ####
  
  
  
  col <- factor(islet_id_interest$Histological_phenotype, levels = unique(islet_id_interest$Histological_phenotype))
  
  
  graphics <- list()
  
  
  
  graphics[[1]] <- function(x, y, w, h){
    grid.rect(x, y, w, h, gp = gpar(fill =  colors[1]))
  }
  
  graphics[[2]] <- function(x, y, w, h){
    grid.rect(x, y, w, h, gp = gpar(fill =  colors[2]))
  }
  
  
  #######
  names(graphics) <- unique(islet_id$Histological_phenotype)
  
  
  
  
  
  
  genesymbol_hm <- expression_list_donor_cat$genesymbol
  
  ## this is the final heatmap
  categories = expression_list_donor_cat$Category
  
  
  if (C==T){
    
    dend = cluster_within_group(as.matrix(x), islet_id_interest$Clinical_phenotype)
    C <- dend
    split <- length(levels(islet_id_interest$Clinical_phenotype))
  }else{
    
    split <- islet_id_interest$Clinical_phenotype
  }
  
  
  
  
  col_fun = colorRamp2(c(min, 0, max), c("blue", "white", "red"), space = "RGB")
  
  if(by_donor == TRUE){
    
    
    
    # all islets
    # all islets
    hm <- Heatmap(as.matrix(x), cluster_rows = R, cluster_columns = C, show_column_names = T, show_row_names =T,
                  ,row_labels = genesymbol_hm, column_names_side = "top", show_column_dend = F, show_row_dend = F
                  # , top_annotation = HeatmapAnnotation(Histo_phenotype = anno_customize(islet_id_interest$Histological_phenotype, graphics = graphics), 
                  #                                      annotation_name_gp = gpar(fontsize = 9))
                  #,cluster_column_slices = TRUE
                  ,column_dend_reorder =  colMeans(as.matrix(x)[1:4,])
                  ,column_split = split, column_labels = islet_id_interest$ID_number
                  ,row_split =categories
                  ,col=col_fun, heatmap_legend_param = list(title="Z-Score",labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=11,fontface="bold"))
                  , column_title_side = "bottom", column_title_rot=0, column_title_gp = gpar(fontsize=12,fontface = "bold"), row_title_rot = 0
                  , row_title_gp = gpar(fontsize = rowtitle
                                        #, fontface = "bold"
                                        ), column_gap = unit(3, "mm") ,row_gap = unit(2, "mm")
                  , row_names_gp = gpar(fontsize = rowname, fontface = "italic"), border = TRUE
                  , column_names_gp = gpar(fontsize = 11), column_names_rot = 60)
    
    # lgd = Legend(title = "Histology",labels_gp = gpar(fontsize = 11), title_gp = gpar(fontsize = 12, fontface="bold"), at = names(graphics), graphics = graphics)
    
    tiff(paste("heatmap_genes_", group_of_interest,"_","by_donor_", process, ".tiff", sep = ""), units="cm", width=width, height=height, res=res)
    
    
    draw(hm, legend_title_gp = gpar(fontsize=11, fontface="bold"))
    #decorate_annotation("Log2FC", {  just = "bottom"})
    
    dev.off()
    
    
  }else{
    hm <- Heatmap(as.matrix(x), cluster_rows = R, cluster_columns = C, show_column_names = F, show_row_names =T,
                  ,row_labels = genesymbol_hm, column_names_side = "top", show_column_dend = F, show_row_dend = F
                   , top_annotation = HeatmapAnnotation(Histology = anno_customize(islet_id_interest$Histological_phenotype, graphics = graphics), 
                                                       annotation_name_gp = gpar(fontsize = 11, fontface="bold"))
                  #,cluster_column_slices = TRUE
                  ,column_split = split
                  #column_labels = islet_id_interest$ID_number
                  
                 # ,column_dend_reorder =  colMeans(as.matrix(x)[1:4,])
                  ,row_split =categories
                  ,col=col_fun, heatmap_legend_param = list(title="Z-Score",labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=11,fontface="bold"))
                  , column_title_side = "bottom", column_title_rot=0, column_title_gp = gpar(fontsize=12,fontface = "bold"), row_title_rot = 0
                  , row_title_gp = gpar(fontsize = rowtitle
                                        #,fontface = "bold"
                                        ), column_gap = unit(3, "mm") ,row_gap = unit(2, "mm")
                  , row_names_gp = gpar(fontsize = rowname, fontface = "italic"), border = TRUE
                  , column_names_gp = gpar(fontsize = 11), column_names_rot = 60)
    
    lgd = Legend(title = "Histology",labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 11, fontface="bold"), at = names(graphics), graphics = graphics)
    
    tiff(paste("heatmap_genes_", group_of_interest,"_","by_islet_", process, ".tiff", sep = ""), units="cm", width=width, height=height, res=res)
    
    
    draw(hm, legend_title_gp = gpar(fontsize=11, fontface="bold"),annotation_legend_list = lgd )
    #decorate_annotation("Log2FC", {  just = "bottom"})
    
    dev.off()
    
  }
  
  
}
