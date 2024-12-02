# set up your working space (folder containing all files), then load all packages and custom functions

# libraries
library(cmapplot)
library(devtools)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gplots)
library(ComplexHeatmap)
library(reshape2)
library(ggradar)
library(ggbreak) 
library(patchwork)
library(plotrix)
library(RColorBrewer)
library(ggfortify)
library(factoextra)
library(plotly)
library(cluster)
library(randomcoloR)
library(scales)
library(tidyverse) 
library(ggpubr)
library(ggplotify)
library(Hmisc)
library(corrplot)
library(stringr)
library(ggrepel)

library(clusterProfiler)
library('org.Hs.eg.db')
library(enrichplot) 
library(ggupset)
library(GOSemSim)
library(pathview)
library(aPEAR)


library(biomaRt)
library(piano)
library(limma)
library(affy)
library (M3C)
library(vsn)
library(GOexpress)
library(EnhancedVolcano)
library(ggVennDiagram)
library(ggvenn)
library(pheatmap)



## loading edited function for radar plots
source(file = "ggradar_edit.R")
environment(ggradar_edit) <- environment(ggradar)

source("heatplot2.R") #vertical heatmaps
environment(heatplot.enrichResult) <- environment(heatplot)


## gene expression islets functions
source("table_ready.R")
source("tukey_tests.R")
source("gene_expression_graphs.R")
source("double_expression_plot.R")
source("multiple_expression_plot.R")
source("scale_values.R")

## pca plots
source("pca_plots.R")

## pca umap plots with M3 package
source("PCA_UMAP_plots.R")

## gene list UMAP
source("gene_list_UMAP.R")

## functions to change labels to words/symbols
source("change_to_words.R")
source("change_to_signs.R")
source("change_to_groups.R")

# to get the DEG results table
source("DEG_data_table.R")
source("barplot_DEGs.R")

# common DEGs
source("common_DEGs.R")

# to make custom volcano plots
source("volcano_plots.R")

# to make custom venn plots
source("venn_plots.R")

# GSEA analysis
source("GO_KEGG_gsea.R")

# functions to visualize GSEA results
source("lollipop_plot.R")
source("pathway_network.R")
source("barplot_result.R")
source("lollipop_plot_enrich.R")
source("barplot_nes.R")

## getting genes from enriched pathways
source("genes_in_enrich.R")

### human to mouse analysis
source("human_to_mouse_genes.R")
source("matching_mouse_human.R")

