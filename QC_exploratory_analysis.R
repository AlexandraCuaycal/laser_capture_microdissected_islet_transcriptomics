########### For exploratory analysis  ############

# create expression set
my_ExpressionSet <- ExpressionSet(assayData=exprs, phenoData = metdata)
save(my_ExpressionSet, file = "my_ExpressionSet.RData")

## quick QC to check data normalization with variance plot
meanSdPlot(exprs(my_ExpressionSet), ranks = TRUE) #data is skewed to the right

## VSN transformation to stabilize variance
glog_data <- justvsn(my_ExpressionSet)
save(glog_data, file = "glog_data.RData")

# VSN transformed data
meanSdPlot(exprs(glog_data), ranks = TRUE) # variance is stabilized

# subset for coding probes only
glog_data_coding <- glog_data[rownames(glog_data) %in% rownames(gene_description_coding), ]
save(glog_data_coding, file = "glog_data_coding.RData")

######################## Making PCA and UMAP plots ###########################

##### PCA plots with untransformed data #####

# variables for plot
cat_var <- varLabels(my_ExpressionSet); cat_var
cat_var <- cat_var[c(1,2,3,4,9)]; cat_var

## saving in new folder "PCA untransformed data"
# use custom function

#can add ellipses with frame=TRUE
pca_plots(my_ExpressionSet, cat_var = cat_var, shape_var = "Clinical_phenotype",folder_name = "PCA untransformed data")

#### PCA plots with transformed data #####
#can add ellipses with frame=TRUE
pca_plots(glog_data, cat_var = cat_var, shape_var = "Clinical_phenotype",folder_name = "PCA VSN transformed data")

###### PLots with M3C package (PCA and UMAP) ############
# uses custom function
## untransformed data
PCA_UMAP_plots(my_ExpressionSet, cat_var = cat_var, folder_name = "Plots untransformed data with M3C package/260 islets")

# VSN transformed data
PCA_UMAP_plots(glog_data, cat_var = cat_var, folder_name = "Plots VSN transformed data with M3C package/260 islets")

#For Huber et al. paper (different color palette)
PCA_UMAP_plots(glog_data, cat_var = cat_var, folder_name = "Plots VSN transformed data with M3C package/260 islets/Huber paper", seed = 45,
               axistextsize = 7, legendtextsize = 7,width = 5, height = 4)

## plots with coding genes only 
# result is similar plot as with glog_data
PCA_UMAP_plots(glog_data_coding, cat_var = cat_var, folder_name = "Plots VSN transformed data with M3C package/260 islets/Coding")


################### removing islets that are enriched in alpha cells ###############

# Run heatmap for cell specific genes first before continuing here

islets <- colnames(glog_data_coding)

## islets with alpha cell enriched probes
load(file = "donors_contam_all_islets.RData")
donors_contam <- donors_contam_all_islets

##islets to remove
islets_rm <- islets[islets %in% names(donors_contam)]; islets_rm

# subsetting variance-stabilized dataset
glog_data_subset<- glog_data[, !(colnames(glog_data) %in% islets_rm)]
save(glog_data_subset, file = "glog_data_subset.RData")

## keeping coding genes only
glog_data_subset_coding <- glog_data_subset[rownames(glog_data_subset) %in% rownames(gene_description_coding), ]
save(glog_data_subset_coding, file = "glog_data_subset_coding.RData")

## Subsetting for islets with HLA Class II risk/neutral
glog_data_subset_HLA_risk_neut<- glog_data_subset[, glog_data_subset$HLA_category %in% c("Risk", "Neutral")]
save(glog_data_subset_HLA_risk_neut, file = "glog_data_subset_HLA_risk_neut.RData")

### subsetting for coding genes only
## this is to be used in DE analysis with piano package
glog_data_subset_HLA_coding<- glog_data_subset_HLA_risk_neut[rownames(glog_data_subset_HLA_risk_neut) %in% rownames(gene_description_coding),]
save(glog_data_subset_HLA_coding,file = "glog_data_subset_HLA_coding.RData")




