install.packages("devtools", dependencies=T)
install.packages("BiocManager", dependencies=T)
install.packages("tidyverse", dependencies=T)
install.packages('Seurat', dependencies=T)
BiocManager::install(c("SingleR","monocle", "DESeq2"),ask = F,update = F) 
BiocManager::install(c("clusterProfiler","DOSE","pheatmap"),ask = F,update = F)
BiocManager::install(c("org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"),ask = F,
                     update = F)

devtools::install_github('RGLab/MAST', upgrade=F, build_vignettes = F)
library(Seurat)
library(MAST)


BiocManager::install(c("SingleCellExperiment"),ask = F,
                     update = F)
