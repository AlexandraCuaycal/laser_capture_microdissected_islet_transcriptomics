###################################################################
#
##### creating heatmap with marker genes for islet cell types
#
#####################################################################

# download gene markers_cell_type_extended.csv from "Other files" folder first

# load islet metdata
islet_id <- read.csv("Islets_Metdata.csv", row.names=1, header=TRUE)
islet_id$Clinical_phenotype <- factor(islet_id$Clinical_phenotype, levels = c("ND", "sAAb", "mAAb", "T1D"))
islet_id$Histological_phenotype <- factor(islet_id$Histological_phenotype, levels = c("INS+CD3-", "INS+CD3+", "INS-CD3-", "INS-CD3+"))

islet_id$Islet_N <- rownames(islet_id)

# rearranging rows
islet_id<- islet_id %>% group_by(Clinical_phenotype) %>%
  arrange(Histological_phenotype)

islet_id <- islet_id[order(islet_id$Clinical_phenotype),]

rownames(islet_id) <- islet_id$Islet_N

# VSN transformed data
load("glog_data_coding.RData")

# expression data 
## all islets 
gene_expression <- exprs(glog_data_coding)
gene_expression<- gene_expression[, match(islet_id$Islet_N, colnames(gene_expression))]

## ordering donors
Donor_ID_levels<- factor(islet_id$Donor_ID, levels = unique(islet_id$Donor_ID))
islet_id<- islet_id[order(unlist(sapply(islet_id$Donor_ID, function(x) which(levels(Donor_ID_levels) == x)))),]

islet_id$ID_number <- factor(islet_id$ID_number, levels = unique(islet_id$ID_number))

### load list of genes for heatmap
gene_list <- read.csv("gene markers_cell_type_extended.csv")

## factor variables
cell_types <- factor(gene_list$Cell_type, levels = c("Beta", "Alpha", "Delta", "Gamma+Epsilon","Ductal", 
                                                      "Acinar", "Stellate", "Endothelial", "Mast", "Macrophage"))

genes <- factor(gene_list$genesymbol, unique(gene_list$genesymbol))

# adding to gene_list
gene_list$Cell_type <- cell_types
gene_list$genesymbol <- genes

                                        
## getting probe_id from gene_description 
## load gene description data
gene_description_data <- read.csv ("Gene_description_data_piano_updated.csv", row.names = 1)

# coding probes only
gene_description_coding <- gene_description_data[gene_description_data$locus.type == "Coding",]
gene_description_coding <- gene_description_coding[grepl("ENSG", gene_description_coding$ensembl_gene_id),]

# subset for cell specific genes
description_list <- gene_description_coding[trimws(gene_description_coding$genesymbol) %in% gene_list$genesymbol,]

# create new data frame                 
probe_list <- data.frame(genesymbol=description_list$genesymbol, probe_id = rownames(description_list))

##ordering probe_id based on gene_list
ordered<- NULL
categories <- cell_types
for (i in 1:length(levels(categories))) {
  
  cat <- levels(categories)
  genes_cat <- gene_list[gene_list$Cell_type == cat[i],]
  gn <- NULL
  df<- NULL
  cat2<- NULL
  
  for (j in 1:length(genes_cat$genesymbol)) {
    g <- probe_list[trimws(probe_list$genesymbol) == genes_cat$genesymbol[j],]
    gn <- rbind(gn, g)
  }
  
  cat2 <- rep(cat[i], times = length(gn$genesymbol))
  df <- data.frame(gn, Category = cat2)
  ordered <- rbind(ordered, df)
}

rownames(ordered) <- ordered$probe_id

# subsetting gene expression with cell specific gene probes
expression_list <- gene_expression[rownames(gene_expression) %in% rownames(ordered),]
expression_list<- expression_list[match(rownames(ordered), rownames(expression_list)),]

# create new dataframe
expression_list_cat <- data.frame(expression_list, genesymbol=ordered$genesymbol, Category=ordered$Category)

# get the genesymbol vector for the heatmap
genesymbol_hm <- expression_list_cat$genesymbol

################### HEATMAP ####################################
library(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)

## for heatmap:
categories <- factor(ordered$Category, levels = unique(ordered$Category)) #this is cell types
genesymbol_hm <- ordered$genesymbol #cell specific genes

# to calculate z scores
x <- expression_list_cat[,1:260] # get only expression values

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


#### for coloring columns based on histo phenotype
histology <- islet_id$Histological_phenotype
                                        
### colors for HISTO phenotype
graphics2 <- list()

graphics2[[1]] <- function(x, y, w, h){
  grid.rect(x, y, w, h, gp = gpar(fill =  color_pal[1]))
}

graphics2[[2]] <- function(x, y, w, h){
  grid.rect(x, y, w, h, gp = gpar(fill =  color_pal[2]))
}

graphics2[[3]] <- function(x, y, w, h){
  grid.rect(x, y, w, h, gp = gpar(fill =  color_pal[3]))
}

graphics2[[4]] <- function(x, y, w, h){
  grid.rect(x, y, w, h, gp = gpar(fill =  color_pal[4]))
}


names(graphics2) <- unique(islet_id$Histological_phenotype)


# to group islets based on Clinical phenotype

col <- factor(islet_id$Clinical_phenotype, levels = unique(islet_id$Clinical_phenotype))

# heatmap z score colors
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "RGB")

## Create heatmap with Complex heatmap
## histo phenotype
hm <- Heatmap(as.matrix(x), cluster_rows = F, cluster_columns = F, show_column_names = F
              ,row_labels = genesymbol_hm, column_names_side = "top"
              , top_annotation = HeatmapAnnotation(Histology = anno_customize(histology, graphics = graphics2), annotation_name_gp = gpar(fontsize = 10))
              ,column_split = islet_id$Clinical_phenotype, row_split =categories, col=col_fun, heatmap_legend_param = list(title="Z-Score",labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10,fontface="bold"))
              , column_title_side = "bottom", column_title_rot=0, column_title_gp = gpar(fontsize=10,fontface = "bold"), row_title_rot = 0
              , row_title_gp = gpar(fontsize = 10,fontface = "bold"), column_gap = unit(2, "mm") ,row_gap = unit(1, "mm")
              , row_names_gp = gpar(fontsize = 5, fontface = "italic"), border = TRUE
              , column_names_gp = gpar(fontsize = 5), column_names_rot = 60)

lgd = Legend(title = "Histology",labels_gp = gpar(fontsize = 9), title_gp = gpar(fontsize = 10, fontface="bold"), at = names(graphics2), graphics = graphics2)

draw(hm, legend_title_gp = gpar(fontsize=8, fontface="bold"), annotation_legend_list = lgd)

                                        
#################### Alpha cell contamination removal ####################################
# heatmap of cell specific genes showed some islets with higher expression of alpha cell probes
# we will use z-scores to exclude those

## getting columns with z-score above 1 sd

z_scores<- x

cell_genes<- expression_list_cat[,261:262]

z_scores_cells <- cbind(z_scores, cell_genes)

## getting only other endocrine cells values
endocrine_cells <- subset(z_scores_cells, Category =="Alpha" | Category== "Delta" |Category== "Gamma+Epsilon")

## getting those islets with 3 or more genes with z-score >1
endocrine<- unique(endocrine_cells$Category)

each_cell_contam_logical<- NULL
for(i in 1:length(endocrine)){
  
  each_cell <- subset(endocrine_cells, Category == endocrine[i])
  each_cell_contam <- apply(each_cell[,1:260],2, function(x){
    sum(x>1)
  })
  
  each_cell_contam_logical <- rbind(each_cell_contam_logical, each_cell_contam)
  
}
rownames(each_cell_contam_logical) <- endocrine

## getting donors with alpha cell enriched islets
donors_contam_all_islets <-  which(each_cell_contam_logical[1,] >= 3); length(donors_contam_all_islets)

## which islets have enriched alpha-cell genes from INS+ islets
islets_INS <- islet_id %>% filter(INS_stain == "Pos") %>% select(Islet_N)

islets_INS_alpha <- donors_contam_all_islets[ names(donors_contam_all_islets) %in% islets_INS$Islet_N]; length(islets_INS_alpha)


save(islets_INS_alpha, file = "islets_INS_alpha.RData")




