library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
dir.create("07_visual")
load(file="./2.results/04_subCl_scRNA.rdata") 
# RidgePlot山脊图
p1 = RidgePlot(scRNA, features = "FCN1")
p2 = RidgePlot(scRNA, features = "PC_2")
plotc = p1/p2 + plot_layout(guides = 'collect')
ggsave('07_visual/ridgeplot_eg.png', plotc, width = 8,height = 8)

# VInPlot小提琴图
p1 = VlnPlot(scRNA, features = "nCount_RNA", pt.size = 0)
p2 = VlnPlot(scRNA, features = "CD8A", pt.size = 0)
plotc = p1/p2 + plot_layout(guides = 'collect')
ggsave('07_visual/vlnplot_eg.png', plotc, width = 8,height = 8)

# FeaturePlot特征图
p1 <- FeaturePlot(scRNA,features = "CD8A", reduction = 'umap')
p2 <- FeaturePlot(scRNA,features = "CD79A", reduction = 'umap')
plotc = p1|p2
ggsave('07_visual/featureplot_eg.png', plotc, width = 10, height = 4)

# DotPlot点图
genelist = c('LYZ','CD79A','CD8A','CD8B','GZMB','FCGR3A')
p = DotPlot(scRNA, features = genelist)
ggsave('07_visual/dotplot_eg.png', p, width = 7, height = 5)

# DoHeatmap热图
genelist = read.csv("03_cellIdentity/top10_diff_genes_wilcox.csv")
genelist <- pull(genelist, gene) %>% as.character
p = DoHeatmap(scRNA, features = genelist, group.by = "seurat_clusters")
ggsave('07_visual/heatmap_eg.png', p, width = 12, height = 9)

# FeatureScatter散点图
p1 <- FeatureScatter(scRNA, feature1 = 'PC_1', feature2 = 'PC_2')
p2 <- FeatureScatter(scRNA, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
plotc = p1|p2
ggsave('07_visual/featurescatter_eg.png', plotc, width = 10, height = 4)

# DimPlot降维图
p1 <- DimPlot(scRNA, reduction = 'tsne', group.by = "celltype", label=T)
p2 <- DimPlot(scRNA, reduction = 'umap', group.by = "Phase", label=T)
p3 <- DimPlot(scRNA, reduction = 'pca', group.by = "celltype", label=T)
p4 <- DimPlot(scRNA, reduction = 'umap', group.by = "seurat_clusters", label=T)
plotc = (p1|p2)/(p3|p4)
ggsave('07_visual/dimplot_eg.png', plotc, width = 10, height = 8)