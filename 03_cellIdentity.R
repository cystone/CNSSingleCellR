library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
rm(list=ls())
dir.create("03_cellIdentity")
load('./2.results/02_cluster_scRNA.rdata')

##以下方法三选一，建议第一种
#默认wilcox方法
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "03_cellIdentity/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "03_cellIdentity/top10_diff_genes_wilcox.csv", row.names = F)
# #专为单细胞设计的MAST
# diff.mast = FindAllMarkers(scRNA, test.use = 'MAST')
# all.markers = diff.mast %>% select(gene, everything()) %>% subset(p_val<0.05)
# top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(all.markers, "03_cellIdentity/diff_genes_mast.csv", row.names = F)
# write.csv(top10, "03_cellIdentity/top10_diff_genes_mast.csv", row.names = F)
# #bulkRNA经典方法DESeq2
# diff.deseq2 = FindAllMarkers(scRNA, test.use = 'DESeq2', slot = 'counts')
# all.markers = diff.deseq2 %>% select(gene, everything()) %>% subset(p_val<0.05)
# top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(all.markers, "03_cellIdentity/diff_genes_deseq2.csv", row.names = F)
# write.csv(top10, "03_cellIdentity/top10_diff_genes_deseq2.csv", row.names = F)
##top10基因绘制热图
top10_genes <- read.csv("03_cellIdentity/top10_diff_genes_wilcox.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave("03_cellIdentity/top10_markers.pdf", plot=plot1, width=8, height=6) 
ggsave("03_cellIdentity/top10_markers.png", plot=plot1, width=8, height=6)


#挑选部分基因
select_genes <- c('LYZ','CD79A','CD8A','CD8B','GZMB','FCGR3A')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0,  ncol=2)
# p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=2)
ggsave("03_cellIdentity/selectgenes_VlnPlot.png", p1, width=6 ,height=8)
#featureplot展示
p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label=T, ncol=2)
ggsave("03_cellIdentity/selectgenes_FeaturePlot.png", p2, width=8 ,height=12)
p3=p1|p2
ggsave("03_cellIdentity/selectgenes.png", p3, width=10 ,height=8)


#----------SingleR鉴定细胞类型----------
library(SingleR)
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"03_cellIdentity/celltype_singleR.csv",row.names = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  ind = which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i])
  scRNA@meta.data[ind,'celltype'] <- celltype$celltype[i]}


p1 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("03_cellIdentity/tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("03_cellIdentity/UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("03_cellIdentity/celltype.pdf", p3, width=10 ,height=5)
ggsave("03_cellIdentity/celltype.png", p3, width=10 ,height=5)

save(scRNA, file="./2.results/03_cell_id_scRNA.rdata")
