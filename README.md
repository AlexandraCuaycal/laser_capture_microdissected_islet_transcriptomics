# Laser-capture microdissected islet transcriptomics
This repository contains scripts and files related to gene expression analyses of laser-capture microdissected islets from nPOD donors across four clinical phenotypes: non diabetic (ND), single auto-antibody positive (sAAb), multiple auto-antibody positive (mAAb) and Type 1 diabetes (T1D) patients. Fresh frozen pancreas slides were immunostained for insulin (INS), CD3, HLA class I and ATP5B to delineate islets with residual beta cells, T cell infiltration, antigen presentation and mitochondrial function. Gene expression analyses were performed with INS+CD3- and INS+CD3+ islets.

## Analysis pipeline
1: Set up working space by loading packages and custom functions.

2: Load expression matrix and metdata.

3: Perform QC and exploratory analysis.

4: Heatmap of islet cell-specific genes.

5: Continue QC (remove alpha cell-enriched islets based on previous step).

6: Get gene expression matrix for genes of interest (Huber et al. paper).

7: Differential expression (DE) analysis.

8: Gene set enrichment analysis (GSEA) for Gene ontology and KEGG.

9: Plot results from GSEA.

10: Gene expression plots for DE(G) genes of interest. 
