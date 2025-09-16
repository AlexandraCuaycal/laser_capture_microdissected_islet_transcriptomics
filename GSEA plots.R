## load cluster profiler GSEA files

load("sAAb_vs_ND_NI_enrichment.RData")
load("T1D_vs_mAAb_NI_enrichment.RData")

load("mAAb_vs_ND_NI_enrichment.RData")
load("T1D_vs_ND_NI_enrichment.RData")

load("unique_sAAb_NI_enrichment.RData")
load("common_sAAb_mAAb_NI_enrichment.RData")
load("unique_mAAb_NI_enrichment.RData")
load("unique_T1D_NI_enrichment.RData")

load("common_T1D_NI_I_enrichment.RData")
load("T1D_vs_mAAb_I_enrichment.RData")
load("T1D_vs_ND_I_vsND_enrichment.RData")
load("T1D_vs_ND_I_upreg_enrichment.RData")

## generating different type of plots
## use the code with the different cluster profiler objects

## for example, using common_sAAb_mAAb_NI_enrichment:

## generating pathway network with GO BP enrichment
pathway_network(common_sAAb_mAAb_NI_enrichment$Readout_tables$GO_enrich@result, title="Shared DEGs sAAb & mAAb vs ND",
                enrichment = T, showcat = 100, gene_number = 10, label_size = 3.5)

## KEGG pathway for Protein processing in endoplasmic reticulum
dme <- pathview(gene.data=common_sAAb_mAAb_NI_enrichment$kegg_gene_list, pathway.id="hsa04141",cpd.idtype = "kegg",
                low = list(gene = "cornflowerblue", cpd = "blue"), mid =
                  list(gene = "gray", cpd = "white"), high = list(gene = "indianred2", cpd =
                                                                    "yellow"), 
                kegg.native=T, same.layer = F, res = 300)

## lollipop plot with GSEA
lollipop_plot(common_sAAb_mAAb_NI_enrichment$GO_KEGG$KEGG_GSEA, title ="Shared DEGs -\nsAAb & mAAb vs ND -\nGSEA KEGG pathways",
              showcat = 10, label_size = 18,cutstring = 25, gene_number = 40 )

## lollipop plot with over-representation result
lollipop_plot_enrich(result =  common_sAAb_mAAb_NI_enrichment$GO_KEGG$GO_enrich, 
                     title = "Shared DEGs -\nsAAb & mAAb vs ND -\nGO BP enrichment",
                     showcat = 9, label_size = 16,cutstring = 40, color = "black" , bg=color_pal[3])

## more plots with GO/KEGG enrichemnts
edox <- setReadable(common_sAAb_mAAb_NI_enrichment$GO_KEGG$GO_enrich, 'org.Hs.eg.db') 

## to make gene network plots with selected categories
edox@result <- edox@result %>% dplyr::filter(Description %in% c(
                                                                  "translational initiation",
                                                                  "positive regulation of translation"
                                                            
                                                                 ))
# need logFC values to add to plot
logFC <- common_sAAb_mAAb_NI_enrichment$logFC

# making the plot
p1 <-cnetplot(edox, cex_label_gene=0.8, cex_label_category=1.2, circular = F, colorEdge = F,
              foldChange = logFC, color_category = "deepskyblue", showCategory =10)

p1 + scale_color_gradient2(low="blue", mid="white", high="red", 
                           limits = c(-2, 2), oob = scales::squish, name = "logFC")+
  scale_size_area(name=str_wrap("Number of genes", width = 10), limits=c(0,100), breaks=c(20,40,60, 80, 100))+
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 11))



