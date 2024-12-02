### function to produce the GO and KEGG GSEA objects and export data
GO_KEGG_gsea <- function(data, folder=NULL, term= "BP",
                         minGO=40, maxGO=150,
                         minKEGG=40, maxKEGG=500){
  
  ## setting organism
  organism <- org.Hs.eg.db
  
  
  df <- data
  
  
  ## continue
  
  df <- df[, c(1,2,3,6,7, 9)]
  
  ## removing duplicated ensembl
  df <- df[!duplicated(df[c("ensembl_gene_ID")]),]
  
  
  # Annotate according to differential expression
  df <- df %>% mutate(diffexpressed = case_when(
    logFC > 0  ~ 'UP',
    logFC < 0  ~ 'DOWN'
  ))
  
  # Entrez ids for KEGG
  
  ensembl_DEGs <- df$ensembl_gene_ID
  
  entrez_DEGs <- bitr(ensembl_DEGs, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
  entrez_DEGs <- entrez_DEGs[!duplicated(entrez_DEGs[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 = df[df$ensembl_gene_ID %in% entrez_DEGs$ENSEMBL,]
  df2 = df2[!duplicated(df2[c("ensembl_gene_ID")]),]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$entrez <- entrez_DEGs$ENTREZID
  
  ## logFC vector
  logFC <-df2$logFC
  names(logFC) <- df2$GeneName
  
  print("ENTREZ mapping done")
  
  ########### GO ###############
  
  gene_list <- df2$logFC
  
  names(gene_list) <- df2$ensembl_gene_ID
  
  kegg_gene_list<- gene_list
  names(kegg_gene_list) <- df2$entrez
  
  ## without logFC
  go_enrich <- enrichGO(gene=df2$ensembl_gene_ID, ont = term, keyType = "ENSEMBL",
                        pvalueCutoff = 0.01, 
                        #qvalueCutoff = 0.01,
                        minGSSize = minGO, #40
                        maxGSSize = maxGO, #200
                        universe = ensembl_gene_id$ensembl_gene_id,
                        OrgDb = organism, 
                        pAdjustMethod = "BH"
  )
  print("GO enrich done")
  
  # with logFC
  go_gsea <- gseGO(geneList = gene_list, ont = term, keyType = "ENSEMBL",
                   pvalueCutoff = 0.01, 
                   minGSSize = minGO, #40
                   maxGSSize = maxGO, #200
                   eps = 0,
                   verbose = TRUE, 
                   OrgDb = organism, 
                   pAdjustMethod = "BH",
                   seed=TRUE
  )
  print("GO GSEA done")
  
  ####### KEGG ANALYISIS #########
  ## kegg gene list
  kegg_gene_list <- df2$logFC
  names(kegg_gene_list) <- df2$entrez
  
  
  ##gene set based on kegg 
  ## considers up and down regulation of genes to identify pathways
  kegg_gsea <- gseKEGG(geneList = kegg_gene_list,
                       pvalueCutoff = 0.01, 
                       minGSSize    = minKEGG, #40
                       maxGSSize    = maxKEGG,
                       keyType       = "ncbi-geneid",
                       eps = 0,
                       verbose = TRUE, 
                       pAdjustMethod = "BH",
                       seed = TRUE
                       
  )
  print("KEGG GSEA done")
  
  ## enrich
  kegg_enrich <- enrichKEGG(gene = df2$entrez,
                            pvalueCutoff = 0.01, 
                            minGSSize    = minKEGG, #40
                            maxGSSize    = maxKEGG,
                            keyType       = "ncbi-geneid",
                            pAdjustMethod = "BH"
  )
  
  print("KEGG GO done")
  
  ### creating list with result objects
  
  cluster_prof_res <- list(GO_enrich = go_enrich, GO_GSEA=go_gsea, KEGG_enrich=kegg_enrich, KEGG_GSEA= kegg_gsea)
  
  ## creating readout tables with genesymbol
  
  readout_tables <- list()
  
  readout_tables_res <- list()
  for(i in 1:length(cluster_prof_res)){
    
    if(grepl( "KEGG",names(cluster_prof_res)[i]) ){
      readout_tables[[i]] <- setReadable(cluster_prof_res[[i]],  'org.Hs.eg.db', 'ENTREZID')
      readout_tables_res[[i]]<- readout_tables[[i]]@result
    }else {
      readout_tables[[i]] <- setReadable(cluster_prof_res[[i]],  'org.Hs.eg.db')
      readout_tables_res[[i]]<- readout_tables[[i]]@result
    }
    
  }
  
  names(readout_tables)<- names(cluster_prof_res)
  names(readout_tables_res) <- names(cluster_prof_res)
  
  print("readout tables done")
  
  if(! is.null(folder) ){
    
    dir.create(folder)
    
    for(i in 1:length(readout_tables_res)){
      write.csv(readout_tables_res[[i]], file = paste(folder, "/",names(readout_tables_res)[[i]], ".csv", sep = ""), row.names = F)
    }
    
  }
  
  return(list(logFC=logFC, gene_list = gene_list, kegg_gene_list=kegg_gene_list, GO_KEGG=cluster_prof_res, Readout_tables=readout_tables))
  
}
