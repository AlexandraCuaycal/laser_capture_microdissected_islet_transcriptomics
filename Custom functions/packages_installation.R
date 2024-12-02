
## install packages needed once

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

BiocManager::install("biomaRt")

BiocManager::install("Orthology.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")

install.packages("aPEAR")
install.packages("reshape2")
install.packages("devtools")
install.packages("ggbreak")
install.packages("patchwork")



## Install current version from GitHub
install.packages("remotes")
remotes::install_github("CMAP-REPOS/cmapplot")


# install.packages("devtools")
devtools::install_github("ricardo-bion/ggradar")


BiocManager::install("limma")
BiocManager::install("piano", dependencies=TRUE)
BiocManager::install("affy")
install.packages("ggfortify")
devtools::install_github("kassambara/factoextra")
install.packages("plotly")
install.packages("randomcoloR")
BiocManager::install("M3C")
BiocManager::install("vsn")
BiocManager::install("GOexpress")
BiocManager::install('EnhancedVolcano')
devtools::install_github("gaospecial/ggVennDiagram")
BiocManager::install("org.Hs.eg.db")
# For visualization
install.packages('pheatmap')
install.packages("DOSE")
install.packages("enrichplot")
install.packages("ggupset")
