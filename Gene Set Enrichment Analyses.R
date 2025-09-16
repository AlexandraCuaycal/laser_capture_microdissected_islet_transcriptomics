######################## GSEA ANALYSIS ############################

## getting gene ontology info for genes in expression data set 

# installing package
install.packages("biomartr")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")

#load packages
library(clusterProfiler)
library(enrichplot)
library(biomartr)
library(biomaRt)

#load objects
load( "ensembl_gene_id_RData")

# DEG objects
load("non_insulitic_data.RData")
load("non_insulitic_data_ND.RData")
load("insulitic_data.RData")
load("insulitic_data_ND.RData")
load("insulitic_vs_non_insulitic_data.RData")

#### start cluster profiler GSEA

################################################################################################
################################# Pseudotimes - not insulitic #################################
#################################################################################################

########################################
############## sAAb vs ND ##############
#########################################
sAAb_vs_ND_NI_enrichment <- GO_KEGG_gsea(data = non_insulitic_data_ND$`sAAb_INS+CD3- vs ND_INS+CD3-`, folder = "sAAb vs ND not insulitic")

save(sAAb_vs_ND_NI_enrichment, file = "sAAb_vs_ND_NI_enrichment.RData")

#######################################################
##################### T1D vs mAAb #####################
######################################################
df<- non_insulitic_data$`T1D_INS+CD3- vs mAAb_INS+CD3-`

T1D_vs_mAAb_NI_enrichment <- GO_KEGG_gsea(data = df, folder = "T1D vs mAAb not insulitic")

save(T1D_vs_mAAb_NI_enrichment, file = "T1D_vs_mAAb_NI_enrichment.RData")


##################################################################################################
########################### Non insulitic - compared to ND ########################################
##################################################################################################

##########################################################
##################### mAAb ##############################
#########################################################
df<- non_insulitic_data_ND$`mAAb_INS+CD3- vs ND_INS+CD3-`

mAAb_vs_ND_NI_enrichment <- GO_KEGG_gsea(data = df, folder = "mAAb vs ND not insulitic")
save(mAAb_vs_ND_NI_enrichment, file = "mAAb_vs_ND_NI_enrichment.RData")


################################################################33333
########################### T1D  ####################################3
####################################################################
df<- non_insulitic_data_ND$`T1D_INS+CD3- vs ND_INS+CD3-`

T1D_vs_ND_NI_enrichment <- GO_KEGG_gsea(data = df, folder = "T1D vs ND not insulitic")
save(T1D_vs_ND_NI_enrichment, file = "T1D_vs_ND_NI_enrichment.RData")

######################################################################################
######################### common genes sAAb and mAAb #################################
######################################################################################

#####################################################
##############  unique in sAAb among all 3 ##########
#####################################################
unique_sAAb_NI_enrichment <- GO_KEGG_gsea(data = unique_sAAb_DEGs,
                                          folder = "unique DEGs in sAAb vs ND NI", )

save(unique_sAAb_NI_enrichment, file = "unique_sAAb_NI_enrichment.RData")

#####################################################
########### genes in sAAb and mAAb ##################
######################################################
# shared genes
common_sAAb_mAAb_NI_enrichment <- GO_KEGG_gsea(data = sAAb_mAAb_DEGs$DEGs_2 %>% filter(Donor_type_expression == "sAAb and mAAb"),
                                          folder = "common DEGs in sAAb and mAAb vs ND NI - done with mAAB logFC")

save(common_sAAb_mAAb_NI_enrichment, file = "common_sAAb_mAAb_NI_enrichment.RData")

#unique genes sAAb and mAAb
unique_sAAb_mAAb_NI_enrichment <- GO_KEGG_gsea(data = unique_sAAb_mAAb_DEGs,
                                               folder = "Unique DEGs in sAAb and mAAb vs ND NI - done with mAAB logFC")

save(unique_sAAb_mAAb_NI_enrichment, file = "unique_sAAb_mAAb_NI_enrichment.RData")

###############################################
########## genes unique mAAb ##################
##############################################
unique_mAAb_NI_enrichment <- GO_KEGG_gsea(data = unique_mAAb_DEGs,
                                               folder = "Unique DEGs mAAb vs ND NI")

save(unique_mAAb_NI_enrichment, file = "unique_mAAb_NI_enrichment.RData")

####################################
########## unique in T1D ############
###################################
unique_T1D_NI_enrichment <- GO_KEGG_gsea(data = unique_T1D_DEGs,
                                          folder = "Unique DEGs T1D vs ND - NI")

save(unique_T1D_NI_enrichment, file = "unique_T1D_NI_enrichment.RData")

##################################################################################################
########################### insulitic contrasts ####################################################
##################################################################################################

####################################################################
##### common genes T1D Not Insulitic and Insulitic ##################
####################################################################

common_T1D_NI_I_enrichment <- GO_KEGG_gsea(data = T1D_NI_I_DEGs$DEGs_1 %>% dplyr::filter(Donor_type_expression == "T1D_NI and T1D_I"),
                                          folder = "Common T1D DEGs NI I vs ND NI")

save(common_T1D_NI_I_enrichment, file = "common_T1D_NI_I_enrichment.RData")

####################################
########## T1D vs mAAb ############
###################################
df<- insulitic_data$`T1D_INS+CD3+ vs mAAb_INS+CD3+`

T1D_vs_mAAb_I_enrichment <- GO_KEGG_gsea(data = df, folder = "T1D vs mAAb insulitic")
save(T1D_vs_mAAb_I_enrichment, file = "T1D_vs_mAAb_I_enrichment.RData")

##################################################################################################
########################### insulitic contrasts - compared to ND #################################
##################################################################################################

####################################
########## T1D vs ND ############
###################################
df<- insulitic_data_ND$`T1D_INS+CD3+ vs ND_INS+CD3-`

T1D_vs_ND_I_vsND_enrichment <- GO_KEGG_gsea(data = df, folder = "T1D insulitic vs ND")
save(T1D_vs_ND_I_vsND_enrichment, file = "T1D_vs_ND_I_vsND_enrichment.RData")

############################################
## getting only upregulated genes in T1D ###
############################################

T1D_insulitic_up <- insulitic_data_ND$`T1D_INS+CD3+ vs ND_INS+CD3-` %>% mutate(diffexpressed = case_when(
  logFC > 0  ~ 'UP',
  logFC < 0  ~ 'DOWN'
))

T1D_insulitic_up <- T1D_insulitic_up %>% filter(diffexpressed == "UP")

T1D_vs_ND_I_upreg <- GO_KEGG_gsea(data = T1D_insulitic_up, folder = "T1D insulitic upregulated vs ND")
save(T1D_vs_ND_I_upreg, file = "T1D_vs_ND_I_upreg_enrichment.RData")
