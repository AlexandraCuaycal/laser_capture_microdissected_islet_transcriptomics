####### To generate expression plots with genes of interest ##############

## This script is used to obtain gene expression data from genes of interest and generate expression plots per-islet based.

#set your working dir to the folder that cointains all files with setwd()

## Huber et al color palette:
## custom color palette for UMAP plot in Huber et al.
color_pal <- c("#000000", "#ea81d7", "#F71482", "#077E97")

## gene expression data
load("glog_data_subset_HLA_coding.RData")
gene_expression <- exprs(glog_data_subset_HLA_coding)

## gene description data
gene_description_data <- read.csv ("Gene_description_data_piano_updated.csv", row.names = 1)
gene_description_data<- gene_description_data %>% dplyr::filter(locus.type == "Coding")
gene_description_data <- gene_description_data[grepl("ENSG", gene_description_data$ensembl_gene_id),]

## islet ID
islet_id <- pData(phenoData(glog_data_subset_HLA_coding))
islet_id$Sample <- rownames(islet_id)

islet_id <- islet_id %>% arrange(Clinical_phenotype)

## to get data only from INS+ CD3- islets
INS_pos_CD3_neg <- islet_id[islet_id$Histological_phenotype == "INS+CD3-",]

#subset gene expression for INS+ CD3- islets
gene_expression_subset_NI <- gene_expression[,colnames(gene_expression) %in% INS_pos_CD3_neg$Sample]

# to get data only from INS + CD3+
INS_pos_CD3_pos <- islet_id[islet_id$Histological_phenotype == "INS+CD3+",]

## removing 1 INS + CD3+ islet from ND donors
INS_pos_CD3_pos <- subset(INS_pos_CD3_pos, Clinical_phenotype != "ND")

#subset gene expression for INS+ CD3+ islets
gene_expression_subset_I <- gene_expression[,colnames(gene_expression) %in% INS_pos_CD3_pos$Sample]


## for convenience, rewrite gene_expression with gene_expression_subset

#Not insulitic (CD3-)
gene_expression_NI <- gene_expression_subset_NI

#Insulitic (CD3+)
gene_expression_I <- gene_expression_subset_I

################ Load genes of interest #########################
## Huber et al. marker genes
genes <- read.csv("Huber et al - target genes.csv")

# get gene description data for genes of interest
genes_list <- gene_description_data[trimws(gene_description_data$genesymbol) %in% genes$genesymbol,]
genes_list <- subset(genes_list, locus.type == "Coding")

## custom function to calculate average, sd and se of expression values
## use expression_values = T to get expression values for all islets

# Get avg expression for INS+CD3- islets
list_expression_NI <- table_ready(gene_expression_NI, genes_list)
write.csv(list_expression_NI, file = "Gene expression avg.csv", row.names = TRUE)

# anova followed by tukey test (islet-based comparisons)
tukey_list_NI <- tukey_tests(gene_expression_NI, genes_list)$tukey_list
write.csv(Tukey_list_NI, file="Tukey tests - genes of interest.csv", row.names = F)


# Create gene expression plots for each gene of interest (islet-based comparisons)
# uses custom function and adds statistical significance based on anova follow by tukey-test (plots comparisons to ND group only)
gene_expression_graphs(gene_expression = gene_expression_NI, 
                       gene_description = gene_description_data, 
                       genes_interest = genes)

############# Donor average comparisons ######################
# Get expression matrix for genes of interest
genes_expression_values <- table_ready(gene_expression_NI, genes_list, expression_values = T)[[1]]
genes_expression_values <- as.data.frame(genes_expression_values)
genes_expression_values$probe_id <- rownames(genes_expression_values)

# Create long expression table
Avg_expression_donors <- pivot_longer(genes_expression_values, values_to = "Expression", names_to = "Sample", cols = 1:78)

# add donor metdata
Avg_expression_donors <- merge(Avg_expression_donors, INS_pos_CD3_neg[,c(1,2,5,6,21)])

#Calculate expression avgs per donor
Avg_expression_donors <- Avg_expression_donors %>% group_by(Donor_ID, probe_id) %>% mutate(Avg_expression_donor = mean(Expression))

# Remove donor avg duplicates
Avg_expression_donors_unique <-subset(Avg_expression_donors, select = - Sample)
Avg_expression_donors_unique <- Avg_expression_donors_unique[! duplicated(Avg_expression_donors_unique[, c("probe_id", "Donor_ID")]),]

# save donor avg expression to csv
write.csv(Avg_expression_donors_unique, file = "Gene expression avg per donor.csv", row.names = F)






