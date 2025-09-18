## create custom heatmaps/polar plots for Genes of Interest (GOIs)
## using custom functions that utilize ComplexHeatmap, ggradar, ggpubr packages

## color palette
color_pal <- c("#000000", "#ea81d7", "#F71482", "#077E97")

# load expression object
load("glog_data_subset_HLA_coding.RData")
gene_expression <- exprs(glog_data_subset_HLA_coding)

## gene description data
gene_description_data <- read.csv ("Gene_description_data_piano_updated.csv", row.names = 1)
gene_description_data<- gene_description_data %>% dplyr::filter(locus.type == "Coding")
gene_description_data <- gene_description_data[grepl("ENSG", gene_description_data$ensembl_gene_id),]

## islet ID
islet_id <- pData(phenoData(glog_data_subset_HLA_coding))
islet_id$Sample <- rownames(islet_id)


## subset islet_id and gene_expression based on islet phenotype

## INS+ CD3- islets
INS_pos_CD3_neg <- islet_id[islet_id$Histological_phenotype == "INS+CD3-",]

### Without T1D islets
INS_pos_CD3_neg_woT1D <- INS_pos_CD3_neg %>% filter(Clinical_phenotype != "T1D")

# INS + CD3+
INS_pos_CD3_pos <- islet_id[islet_id$Histological_phenotype == "INS+CD3+",]
## removing 1 islet from ND donors
INS_pos_CD3_pos <- subset(INS_pos_CD3_pos, Clinical_phenotype != "ND")


gene_expression_subset_NI <- gene_expression[,colnames(gene_expression) %in% INS_pos_CD3_neg$Sample]
gene_expression_subset_I <- gene_expression[,colnames(gene_expression) %in% INS_pos_CD3_pos$Sample]

gene_expression_subset_NI_woT1D <- gene_expression[,colnames(gene_expression) %in% INS_pos_CD3_neg_woT1D$Sample]


## for convenience, rename gene_expression_subset with gene_expression

#Not insulitic
gene_expression_NI <- gene_expression_subset_NI; rm(gene_expression_subset_NI)

#Insulitic
gene_expression_I <- gene_expression_subset_I ; rm(gene_expression_subset_I)


## Insulitic with ND data
islets_ND <- INS_pos_CD3_neg %>% filter(Clinical_phenotype == "ND", Histological_phenotype == "INS+CD3-")
ND_expression <- gene_expression_NI[,colnames(gene_expression_NI) %in% islets_ND$Sample]

insulitic_islets_ND <- rbind(islets_ND, INS_pos_CD3_pos)
insulitic_expression_ND <- cbind(ND_expression, gene_expression_I)

## T1D islets only
T1D_islets <- islet_id %>% filter(Clinical_phenotype == "T1D")

T1D_islets_ND <- rbind(islets_ND, T1D_islets)
T1D_islets_ND$Clinical_phenotype <- factor(T1D_islets_ND$Clinical_phenotype, levels = c("ND", "T1D"))

T1D_expression <- gene_expression[,colnames(gene_expression) %in% T1D_islets$Sample]
T1D_expression_ND <- cbind(ND_expression, T1D_expression)

######################### continuing with example for common_sAAb_mAAb_NI_enrichment: ############################

## to generate heatmap of GOIs ##
## uses custom function to extract genes enriched in specific GO/KEGG and custom heatmap function ##

### ER protein processsing KEGG pathway
ER_sAAb_mAAb<-  genes_in_enrich(enrichment = common_sAAb_mAAb_NI_enrichment$Readout_tables$KEGG_GSEA, 
                                      description = c( 
                                        "Protein processing in endoplasmic reticulum"
                                      ), 
                                      gene_description = gene_description_data,
                                      DEG_table = non_insulitic_data_ND$`mAAb_INS+CD3- vs ND_INS+CD3-`)

ER_protein_transport <- ER_sAAb_mAAb %>% filter(genesymbol %in% 
                                        c("SEC61A", "RRBP1",  "LMAN1", "LMAN2", "SEC31A", "SEC24C"))
ER_protein_transport$Category <- rep("protein transport in/out of ER", times=nrow(ER_protein_transport))

ER_folding <- ER_sAAb_mAAb %>% filter(genesymbol %in% 
                                        c("HYOU1", "P4HB", "CANX", "CALR",
                                          "UGGT1", "PRKCSH"))

ER_folding$Category <- rep("Protein folding/glycosylation in ER", times=nrow(ER_folding))

ER_stress <- ER_sAAb_mAAb %>% filter(genesymbol %in% 
                                       c("VCP", "HSPA6", "ATF6", "ATF6B", "ERN1", "XBP1"))

ER_stress$Category <- rep("ER stress response", times=nrow(ER_stress))

ER_genes <- bind_rows(ER_protein_transport, ER_folding, ER_stress)


## BPs related to translation
Translation_sAAb_mAAb <-genes_in_enrich(enrichment = common_sAAb_mAAb_NI_enrichment$Readout_tables$GO_enrich, 
                                        description = c( 
                                          "translational initiation",
                                          "positive regulation of translation"
                                        ), 
                                        gene_description = gene_description_data,
                                        DEG_table = non_insulitic_data_ND$`mAAb_INS+CD3- vs ND_INS+CD3-`)

Initiation_factors <- Translation_sAAb_mAAb %>% filter(genesymbol %in% c(
  "EIF2A","EIF3C", "EIF3E", "EIF3I", "EIF3M", "EIF4A1", "EIF4B", "EIF4H", "EIF6"
))
Initiation_factors$Category <- rep("Translation initiation factors", times=nrow(Initiation_factors))

Reg_translation <- Translation_sAAb_mAAb %>% filter(genesymbol %in% c(
  "EIF2AK1", "EIF4EBP2", "AGO2", "MTOR", "DDX3X"
))
Reg_translation$Category <- rep("Regulators of translation", times=nrow(Reg_translation))

## ribosomal subunits
Ribosome_40 <- c("RPS2", "RPS3", "RPS9", "RPS25")
small_ribosome <-  non_insulitic_data_ND$`mAAb_INS+CD3- vs ND_INS+CD3-` %>% filter(GeneName %in% Ribosome_40)
small_ribosome$Category <- rep('Small ribosomal proteins', times=nrow(small_ribosome))
colnames(small_ribosome)[1:2] <- c("probe_id", "genesymbol")

Ribosome_60 <- c("RPL4", "RPL8", "RPL10", "RPL10A", "RPL15", "RPL17", "RPL18A", "RPL28")
large_ribosome <-  non_insulitic_data_ND$`mAAb_INS+CD3- vs ND_INS+CD3-` %>% filter(GeneName %in% Ribosome_60)
large_ribosome$Category <- rep('Large ribosomal proteins', times=nrow(large_ribosome))
colnames(large_ribosome)[1:2] <- c("probe_id", "genesymbol")

translation <- bind_rows(Initiation_factors, Reg_translation, small_ribosome, large_ribosome)

##lysosome KEGG pathway
lysosome_sAAb_mAAb <-genes_in_enrich(enrichment = common_sAAb_mAAb_NI_enrichment$Readout_tables$KEGG_GSEA, 
                                        description = c( 
                                          "Lysosome"
                                        ), 
                                        gene_description = gene_description_data,
                                        DEG_table = non_insulitic_data_ND$`mAAb_INS+CD3- vs ND_INS+CD3-`)

lysosome_sAAb_mAAb <- lysosome_sAAb_mAAb %>% filter(genesymbol %in% c("IGF2R", "CLTC","CTSC",  "GBA1", "GNS", "ACP2", 
                                                                      "ASAH1", "PSAP" ,"LAMP2", "SCARB2", "ATP6V0C"))

lysosome_sAAb_mAAb$Category <- rep("Lysosome", times=nrow(lysosome_sAAb_mAAb))

### bind all tables together for heatmap ##
Shared_sAAb_mAAb <- bind_rows(translation, ER_genes, lysosome_sAAb_mAAb)

Shared_sAAb_mAAb$Category <- str_wrap(Shared_sAAb_mAAb$Category, width = 20)
Shared_sAAb_mAAb$Category <- factor(Shared_sAAb_mAAb$Category, levels = unique(Shared_sAAb_mAAb$Category))

## Generate heatmap with custom function
## will export figure to working dir
hm_genes_of_interest_cat(Shared_sAAb_mAAb[!duplicated(Shared_sAAb_mAAb$genesymbol),],
                         gene_description_data = gene_description_data,
                         gene_expression_obj = glog_data_subset_HLA_coding,
                         group_of_interest = c("sAAb","mAAb"),
                         Insulitis_cat = "Not_insulitic", 
                         process = "shared sAAb mAAb",
                         width = 20,
                         height = 22, by_donor = FALSE)


################ To generate radial plots ##############################
###### uses custom function with ggradar package

### for not insulitic islets:
ER_stress_plot <- multiple_expression_plot(gene_expression = gene_expression_NI, islet_id = INS_pos_CD3_neg, 
                                           gene_description = gene_description_data, 
                                        genes_interest = ER_stress %>%
                                          dplyr::select(probe_id, genesymbol), legend_size=16,
                                        title="INS+CD3- islets", title_size=18, label_size = 4.5)

## for insulitic islets: 
ER_stress_plot_CD3 <- multiple_expression_plot(gene_expression = insulitic_expression_ND, islet_id = insulitic_islets_ND, 
                                           gene_description = gene_description_data, 
                                           genes_interest = ER_stress %>%
                                             dplyr::select(probe_id, genesymbol), legend_size=16,
                                           title="INS+CD3+ islets", title_size=18, label_size = 4.5)

## plotting both together with patchwork package
ER_stress_plot + ER_stress_plot_CD3 + plot_layout(guides = "collect", ncol=1) +
  plot_annotation(title = "ER stress", theme = theme(plot.title = element_text(hjust = 0.3, size = 18, 
                                                                               colour = "darkblue", face = "bold")))


ggsave(filename = "ER stress_islets.tiff",height = 22, width = 18, units = "cm", dpi = 300 )

### Correlation plots ####
# uses custom function with ggpubr package

## ER stress markers in non insulitic islets
### ER stress
XBP1_CD68 <- double_expression_plot_corr(gene_expression = gene_expression_NI, islet_id = INS_pos_CD3_neg, 
                                         gene_description = gene_description_data, genes = c("CD68", "XBP1"),
                                         ylim=c(8, 15), color = color_pal, x_title = "")

XBP1_CD3E <- double_expression_plot_corr(gene_expression = gene_expression_NI, islet_id = INS_pos_CD3_neg, 
                                           gene_description = gene_description_data, genes = c("CD3E", "XBP1"),
                                           ylim=c(8, 15), color = color_pal, x_title = "")

ERN1_CD68 <- double_expression_plot_corr(gene_expression = gene_expression_NI, islet_id = INS_pos_CD3_neg, 
                                         gene_description = gene_description_data, genes = c("CD68", "ERN1"),
                                         ylim=c(5, 13), color = color_pal, x_title = "")

ERN1_CD3E <- double_expression_plot_corr(gene_expression = gene_expression_NI, islet_id = INS_pos_CD3_neg, 
                                         gene_description = gene_description_data, genes = c("CD3E", "ERN1"),
                                         ylim=c(5, 13), color = color_pal, x_title = "")


XBP1_CD68 + XBP1_CD3E + ERN1_CD68 + ERN1_CD3E + plot_layout(guides = "collect", axes = "collect")+
  plot_annotation(title = "INS+CD3- islets",
                  theme = theme(plot.title = element_text(size = 17, face = "bold", colour = "darkblue")))

ggsave(filename = "corr_plot_NI_ER stress_CD68_CD3.tiff", width = 12, height = 12, dpi = 300)

## ER stress markers in insulitic islets

XBP1_CD68_I <- double_expression_plot_corr(gene_expression = gene_expression_I, islet_id = INS_pos_CD3_pos, 
                                           gene_description = gene_description_data, genes = c("CD68", "XBP1"),
                                           ylim=c(8, 15), color = color_pal[2:4], x_title = "")

XBP1_CD3E_I <- double_expression_plot_corr(gene_expression = gene_expression_I, islet_id = INS_pos_CD3_pos, 
                                         gene_description = gene_description_data, genes = c("CD3E", "XBP1"),
                                         ylim=c(8,15), color = color_pal[2:4], x_title = "")


ERN1_CD68_I <- double_expression_plot_corr(gene_expression = gene_expression_I, islet_id = INS_pos_CD3_pos, 
                                           gene_description = gene_description_data, genes = c("CD68", "ERN1"),
                                           ylim=c(5, 13), color = color_pal[2:4], x_title = "")

ERN1_CD3E_I <- double_expression_plot_corr(gene_expression = gene_expression_I, islet_id = INS_pos_CD3_pos, 
                                        gene_description = gene_description_data, genes = c("CD3E", "ERN1"),
                                        ylim=c(5, 13), color = color_pal[2:4], x_title = "")



XBP1_CD68_I + XBP1_CD3E_I + ERN1_CD68_I + ERN1_CD3E_I + plot_layout(guides = "collect", axes = "collect")+
  plot_annotation(title = "INS+CD3+ islets",
                  theme = theme(plot.title = element_text(size = 17, face = "bold", colour = "darkblue")))

ggsave(filename = "corr_plot_NI_ER stress_CD68_CD3pos.tiff", width = 12, height = 12, dpi = 300)


#### Gene expression plots across clinical phenotype ###

## uses custom function to plot gene expression of GOI across islets

gene_expression_graphs_sig(gene_expression=gene_expression_NI, 
                                                 gene_description = gene_description_data, 
                                                 genes_interest = "ERN1", folder=NULL, 
                                                 width=16, height=14,
                                                 x_title="", 
                                                 color= color_pal, 
                                                 islet_id= islet_id, legend_size=16,
                                                 axis_ticks_size=14)


