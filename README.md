# Laser-capture microdissected islet transcriptomics
This repository contains scripts and files related to gene expression analyses of laser-capture microdissected islets from donors with no diabetes (ND), single auto-antibody positive (sAAb), multiple auto-antibody positive (mAAb) and Type 1 diabetes (T1D) patients. Fresh frozen pancreas sections were obtained from the Network for Pancreatic Organ donors with Diabetes (nPOD). Pancreas sections were immunostained for insulin (INS), CD3, HLA class I and ATP5B to delineate islets with residual beta cells, T cell infiltration, antigen presentation and mitochondrial function. Gene expression analyses were performed in R Studio with INS+CD3- and INS+CD3+ islets.

## Data availability
CEL files and expression values for each islet can be found at GEO database with GEO identifier: [GSE284772](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284772). 260 islets expression matrix is downloaded from https://doi.org/10.5281/zenodo.14537115

## Analysis pipeline
1: Set up working space by loading packages and custom functions.

2: Load expression matrix and metdata.

3: Perform QC and exploratory analyses.

4: Heatmap of islet cell-specific genes.

5: Continue QC (remove alpha cell-enriched islets based on previous step).

6: Get donor expression averages for genes of interest (Huber et al. paper).

7: Differential expression (DE) analyses.

8: Gene set enrichment analyses (GSEA) with gene ontology and KEGG.

9: Plot results from GSEA.

10: Gene expression plots for DE(G) genes of interest. 
