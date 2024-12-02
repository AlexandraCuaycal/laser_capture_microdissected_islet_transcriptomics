## Start by setting up color palette for plots

## General color palette for analysis
color_pal <- brewer.pal(6, "Dark2")[c(1,3,4,6)]
color_pal<- color_pal[c(1, 4,2,3)]

# visualize colors
show_col(color_pal)

## Huber et al color palette:
## custom color palette for UMAP plot in Huber et al.
color_pal <- c("#000000", "#ea81d7", "#F71482", "#077E97")

##### Load and clean files ###########

## gene expression data
exprs <- as.matrix(read.csv("gene_expression_groups.csv", header=TRUE,
                            row.names=1,
                            as.is=TRUE))

## met data for the Islets
metdata <-  read.csv("Islets_Metdata.csv", row.names=1, header=TRUE)

## reordering columns in exprs to match order in metdata
exprs <- exprs[, rownames(metdata)]

#factor variables in metdata
metdata$Clinical_phenotype<- factor(metdata$Clinical_phenotype, levels = unique(metdata$Clinical_phenotype))
metdata$Histological_phenotype<- factor(metdata$Histological_phenotype, levels = unique(metdata$Histological_phenotype))
metdata$Donor_ID<- factor(metdata$Donor_ID, levels = unique(metdata$Donor_ID))
metdata$Histological<- factor(metdata$Histological, levels = unique(metdata$Histological))
metdata$ID_number<- factor(metdata$ID_number, levels = unique(metdata$ID_number))
metdata$INS_stain<- factor(metdata$INS_stain, levels = unique(metdata$INS_stain))
metdata$Insulitis<- factor(metdata$Insulitis, levels = unique(metdata$Insulitis))
metdata$Race<- factor(metdata$Race, levels = unique(metdata$Race))
metdata$HLA_category<- factor(metdata$HLA_category, levels = unique(metdata$HLA_category))
metdata$Autoantibody_type<- factor(metdata$Autoantibody_type, levels = unique(metdata$Autoantibody_type))
metdata$COD_ID<- factor(metdata$COD_ID, levels = unique(metdata$COD_ID))

metdata<- subset(metdata, select= - Histological)

## converting metdata to Annotated data frame
variables_metdata <- data.frame(Var_description = colnames(metdata))
metdata <- new("AnnotatedDataFrame", data=metdata, varMetadata=variables_metdata)

## metdata for probes is gene_description_data
## gene description data
gene_description_data <- read.csv ("Gene_description_data_piano_updated.csv", row.names = 1)

gene_description_coding <- gene_description_data[gene_description_data$locus.type == "Coding",]

gene_description_coding <- gene_description_coding[grepl("ENSG", gene_description_coding$ensembl_gene_id),]


########################## curating gene_description data ################

## coding genes
gene_description_data <- gene_description_coding


## getting chromosome number
chromosome <- unlist(lapply(strsplit(gene_description_data$chromosome, split = "chr"), function(x) {
  x<-unlist(x)[2]
}))

## long chromosome locations
for (i in 1:length(chromosome)){
  
  if(grepl("_", chromosome[i] ) == TRUE){
    chromosome[i] <- unlist(strsplit(chromosome[i], split="_"))[1]
  }
  
}

# adding to gene description data table
gene_description_data$chromosome <- chromosome

## assigning numbers to sexual, mito and unkown chr
X_chromosomes <- which(gene_description_data$chromosome == "X") 
Y_chromosomes <- which(gene_description_data$chromosome == "Y")
M_chromosomes <- which(gene_description_data$chromosome == "M")
Unknown_chromosomes <- which(gene_description_data$chromosome == "Un")

gene_description_data$chromosome[X_chromosomes]<- 23
gene_description_data$chromosome[Y_chromosomes]<- 23
gene_description_data$chromosome[M_chromosomes]<- 24
gene_description_data$chromosome[Unknown_chromosomes]<- NA

NAs <- which(is.na(gene_description_data$chromosome))
gene_description_data$chromosome[NAs]<- 25

gene_description_data$chromosome<- as.numeric(gene_description_data$chromosome)

NAs <- which(is.na(gene_description_data$start))
gene_description_data$start[NAs]<- 0

## gene description table for piano package with coding genes
gene_description_data_piano <- gene_description_data[, c(1,2,3)]

#ensembl gene IDs for all coding probes
ensembl_gene_id <- gene_description_data[,c(1,7)]
save(ensembl_gene_id, file = "ensembl_gene_id_RData")



