library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
rm(list=ls())
dir.create("04_subCluster")
load(file="./2.results/03_cell_id_scRNA.rdata")
##提取细胞子集
Cells.sub <- subset(scRNA@meta.data, celltype=="T_cells")
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))

# 提重新降维聚类



# 因为再聚类的细胞之间差异比较小，所以聚类函数FindClusters()控制分辨率的参数建议调高到resolution = 0.9。

##PCA降维
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:10
##细胞聚类
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.9)
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'04_subCluster/cell_cluster.csv',row.names = F)
##非线性降维
#tSNE
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
write.csv(embed_tsne,'04_subCluster/embed_tsne.csv')
plot1 = DimPlot(scRNAsub, reduction = "tsne") 
ggsave("04_subCluster/tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("04_subCluster/tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
write.csv(embed_umap,'04_subCluster/embed_umap.csv') 
plot2 = DimPlot(scRNAsub, reduction = "umap") 
ggsave("04_subCluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("04_subCluster/UMAP.png", plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("04_subCluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("04_subCluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)


# Cluster差异分析
diff.wilcox = FindAllMarkers(scRNAsub)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "04_subCluster/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "04_subCluster/top10_diff_genes_wilcox.csv", row.names = F)

# SingleR细胞鉴定

##细胞类型鉴定
library(SingleR)
refdata <- MonacoImmuneData()
testdata <- GetAssayData(scRNAsub, slot="data")
clusters <- scRNAsub@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"04_subCluster/celltype_singleR.csv",row.names = F)
scRNAsub@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("04_subCluster/tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("04_subCluster/UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("04_subCluster/celltype.pdf", p3, width=10 ,height=5)
ggsave("04_subCluster/celltype.png", p3, width=10 ,height=5)