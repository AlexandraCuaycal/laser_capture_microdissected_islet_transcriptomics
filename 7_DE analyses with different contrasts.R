
## for Piano Package (DE and GSEA analysis)
## gene description for coding genes
load("gene_description_data_piano.RData")

## ensembl gene ids
load( "ensembl_gene_id_RData")

#changing genesymbol column to geneName for piano package
colnames(gene_description_data_piano)[1]<- "geneName"

# alpha cell enriched removed from INS+ islets
##  removed HLA risk neutral
load("glog_data_subset_HLA_coding.RData")

## metdata with 260 minus alpha cell enriched & only HLA risk, neutral
mySetup <- pData(phenoData(glog_data_subset_HLA_coding))
mySetup <- mySetup[, 1:2]

# converting symbol to word for piano package
histo_phenotype_words <- change_to_words(string =as.character(mySetup$Histological_phenotype))

# replacing in metdata table
mySetup$Histological_phenotype <- factor(histo_phenotype_words, levels = unique(histo_phenotype_words))

# create array object
myArrayData <- loadMAdata(dataNorm = exprs(glog_data_subset_HLA_coding), annotation = gene_description_data_piano, setup = mySetup, platform = "NULL")

## getting factors to build contrasts
factors_data <- extractFactors(myArrayData)
unique_factors <- unique(factors_data); unique_factors

#### DEG analysis ####

## creating folders for different comparisons

#Non insulitic between groups across pseudotime
folder1 <- "Non Insulitic islets between groups/Between pseudo times"
dir.create(folder1)

#Non insulitic vs control
folder4 <- "Non Insulitic islets between groups/Compared to ND group"
dir.create(folder4)

## Insulitic betw groups
folder2 <- "Insulitic islets between groups"
dir.create(folder2)

## Insulitic vs Non-insulitic within groups
folder3 <- "Insulitic vs non-insulitic within each group"
dir.create("Insulitic vs non-insulitic within each group")

#### not insulitic across pseudotimes
non_insulitic_contrasts <- c( "sAAb_INS_pos_CD3_neg - ND_INS_pos_CD3_neg", "mAAb_INS_pos_CD3_neg - sAAb_INS_pos_CD3_neg", "T1D_INS_pos_CD3_neg - mAAb_INS_pos_CD3_neg")

DEGs_non_insulitic <- diffExp(myArrayData, contrasts = non_insulitic_contrasts, colors = color_pal,
                              volcanoFC = 1.5, significance = 0.001, save = TRUE)

save(DEGs_non_insulitic, file = "DEGs_non_insulitic.RData")


## not insulitic vs the ND group
non_insulitic_contrasts_ND <- c("sAAb_INS_pos_CD3_neg - ND_INS_pos_CD3_neg", "mAAb_INS_pos_CD3_neg - ND_INS_pos_CD3_neg", "T1D_INS_pos_CD3_neg - ND_INS_pos_CD3_neg")

DEGs_non_insulitic_ND <- diffExp(myArrayData, contrasts = non_insulitic_contrasts_ND, colors = color_pal,
                              volcanoFC = 1.5, significance = 0.001, save = TRUE)

save(DEGs_non_insulitic_ND, file = "DEGs_non_insulitic_vsND.RData")

#### insulitic islets
## comparing INS+ CD3+ islets 
insulitic_contrasts <- c( "mAAb_INS_pos_CD3_pos - sAAb_INS_pos_CD3_pos", "T1D_INS_pos_CD3_pos - mAAb_INS_pos_CD3_pos", "T1D_INS_pos_CD3_pos - sAAb_INS_pos_CD3_pos")

##using pval 0.001 for gsea
DEGs_insulitic <- diffExp(myArrayData, contrasts = insulitic_contrasts, colors = color_pal,
                          volcanoFC = 1, significance = 0.001, save = TRUE, heatmapCutoff = 1e-5)

save(DEGs_insulitic, file = "DEGs_insulitic.RData")

#insulitic compared to ND
insulitic_contrasts_ND <- c( "mAAb_INS_pos_CD3_pos - ND_INS_pos_CD3_neg", "T1D_INS_pos_CD3_pos - ND_INS_pos_CD3_neg")

DEGs_insulitic_ND <- diffExp(myArrayData, contrasts = insulitic_contrasts_ND, colors = color_pal,
                          volcanoFC = 1, significance = 0.001, save = TRUE, heatmapCutoff = 1e-5)


save(DEGs_insulitic_ND, file = "DEGs_insulitic_ND.RData")

### comparing specific groups INS+
# comparing insulitic vs non-insulitic for each group (sAAb, mAAb and T1D)

based_on_insulitis_contrasts <- c( "sAAb_INS_pos_CD3_pos - sAAb_INS_pos_CD3_neg", "mAAb_INS_pos_CD3_pos - mAAb_INS_pos_CD3_neg", "T1D_INS_pos_CD3_pos - T1D_INS_pos_CD3_neg")

DEGs_based_on_insulitis <- diffExp(myArrayData, contrasts = based_on_insulitis_contrasts, colors = color_pal,
                                   volcanoFC = 1, significance = 0.001, save = TRUE, heatmapCutoff = 1e-5)

save(DEGs_based_on_insulitis, file = "DEGs_based_on_insulitis.RData")


### generating DEG tables with adj. p val and logFC cutoffs
### NOT INSULITIC ###
##data tables
non_insulitic_data <- DEG_data_table(DEG_object = DEGs_non_insulitic, folder = folder1, 
                                     adj.p.val = 0.001, logFCcutoff = 0.5,
                                     contrast_to_get = non_insulitic_contrasts)

save(non_insulitic_data, file = "non_insulitic_data.RData")

non_insulitic_data_ND <- DEG_data_table(DEG_object = DEGs_non_insulitic_ND, folder = folder4, logFCcutoff = 0.5,
                                     adj.p.val = 0.001, contrast_to_get = non_insulitic_contrasts_ND)

save(non_insulitic_data_ND, file = "non_insulitic_data_ND.RData")

# volcano & venn plots
volcano_plots(data=DEGs_non_insulitic, DE_contrasts = non_insulitic_contrasts,folder_name = folder1, title = "INS+CD3- islets",
              FCcutoff = 0.5, pCutoff = 0.001, labsize = 3.5,
              width=16, height=16, dpi=600, max_overlaps = 2,
              titleLabSize = 16,
              subtitleLabSize = 14,
              captionLabSize = 14,
              legendLabSize = 14,
              legendIconSize = 3,
              axisLabSize = 14)

venn_plots(data=DEGs_non_insulitic, DE_contrasts = non_insulitic_contrasts,folder_name = folder1, 
           logFCcutoff = 0.5, adj.p.val = 0.001,
           title = "INS+CD3- islets - Pseudotime")

#barplot
barplot_DEGs(non_insulitic_data, 
             Contrast = c("sAAb vs ND", "mAAb vs sAAb", "T1D vs mAAb"),
             Title = "INS+CD3- islets - Pseudotime")

## vs ND
volcano_plots(data=DEGs_non_insulitic_ND, DE_contrasts = non_insulitic_contrasts_ND,folder_name = folder4, 
              title = "INS+CD3- islets", pCutoff = 0.001, labsize = 3.5,
              width=16, height=16, dpi=600, max_overlaps = 2,
              titleLabSize = 16,
              subtitleLabSize = 14,
              captionLabSize = 14,
              legendLabSize = 14,
              legendIconSize = 3,
              axisLabSize = 14) 
venn_plots(data=DEGs_non_insulitic_ND, DE_contrasts = non_insulitic_contrasts_ND,folder_name = folder4, 
           logFCcutoff = 0.5, width = 12, height = 12,
           title = "INS+CD3- islets - Compared to ND")

# barplots
barplot_DEGs(non_insulitic_data_ND, 
             Contrast = c("sAAb vs ND", "mAAb vs ND", "T1D vs ND"),
             Title = "INS+CD3- islets - Compared to ND")

### INSULITIC ###

# pval 0.001
##data tables
insulitic_data <- DEG_data_table(DEG_object = DEGs_insulitic, adj.p.val = 0.001, 
                                 folder = folder2, contrast_to_get = insulitic_contrasts)

save(insulitic_data, file ="insulitic_data.RData" )

## volcano and venn plots
volcano_plots(data=DEGs_insulitic, DE_contrasts = insulitic_contrasts, pCutoff = 0.001, FCcutoff = 0.5 ,folder_name = folder2, title = "INS+CD3+ islets")

venn_plots(data=DEGs_insulitic, DE_contrasts = insulitic_contrasts,folder_name = folder2,
           title = "INS+CD3+ islets - Pseudotime")


## compared to ND
insulitic_data_ND <- DEG_data_table(DEG_object = DEGs_insulitic_ND, adj.p.val = 0.001, 
                                    folder = folder2, contrast_to_get = insulitic_contrasts_ND)

save(insulitic_data_ND, file ="insulitic_data_ND.RData" )

#volcano and venn plots
volcano_plots(data=DEGs_insulitic_ND, DE_contrasts = insulitic_contrasts_ND, pCutoff = 0.001, FCcutoff = 0.5 ,folder_name = paste(folder2, "/compared to ND", sep = ""), title = "INS+CD3+ islets")

venn_plots(data=DEGs_insulitic_ND, DE_contrasts = insulitic_contrasts_ND,folder_name = folder2,adj.p.val = 0.001,
           title = "INS+CD3+ islets - Compared to ND", alpha = 0, cat.names = c("mAAb", "T1D"))

# barplots
barplot_DEGs(insulitic_data_ND, 
             Contrast = c("mAAb vs ND", "T1D vs ND"),
             Title = "INS+CD3+ islets - Compared to ND INS+CD3-")


### comparing insulitic vs non-insulitic for each group (sAAb, mAAb and T1D)
##data tables
insulitic_vs_non_insulitic_data <- DEG_data_table(DEG_object = DEGs_based_on_insulitis, adj.p.val=0.001, folder = folder3, contrast_to_get = based_on_insulitis_contrasts)

save(insulitic_vs_non_insulitic_data, file = "insulitic_vs_non_insulitic_data.RData")

# volcano plots
volcano_plots(data=DEGs_based_on_insulitis, DE_contrasts = based_on_insulitis_contrasts,folder_name = folder3, 
              pCutoff = 0.001, FCcutoff = 0.5,
              title = "Based on insulitis within INS+ islets")

venn_plots(data=DEGs_based_insulitis, DE_contrasts = based_on_insulitis_contrasts[2:3],folder_name = folder3,adj.p.val = 0.001,
           title = "INS+CD3+ vs INS+CD3- islets", cat.names = c("mAAb", "T1D"))

#barplot
barplot_DEGs(insulitic_vs_non_insulitic_data[2:3], 
             Contrast = c("mAAb", "T1D"),
             Title = "INS+CD3+ vs INS+CD3- islets")


