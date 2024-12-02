###################################################################
#
##### creating heatmap with marker genes for islet cell types
#
#####################################################################

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


description_list <- gene_description_coding[trimws(gene_description_coding$genesymbol) %in% gene_list$genesymbol,]

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

expression_list <- gene_expression[rownames(gene_expression) %in% rownames(ordered),]

expression_list<- expression_list[match(rownames(ordered), rownames(expression_list)),]

expression_list_cat <- data.frame(expression_list, genesymbol=ordered$genesymbol, Category=ordered$Category)


genesymbol_hm <- expression_list_cat$genesymbol

