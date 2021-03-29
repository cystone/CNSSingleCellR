library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
rm(list=ls())
dir.create("05_pseudoTime")

#--------数据导入与处理---------
load(file="./2.results/04_subCl_scRNA.rdata") #scRNAsub是上一节保存的T细胞子集seurat对象
data <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
#mycds <- detectGenes(mycds, min_expr = 2)  #很多教程不用

# 提选择代表性基因

##使用clusters差异表达基因
diff.genes <- read.csv('./04_subCluster/diff_genes_wilcox.csv')
diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
##使用seurat选择的高变基因
var.genes <- VariableFeatures(scRNA)
mycds <- setOrderingFilter(mycds, var.genes)
p2 <- plot_ordering_genes(mycds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
##结果对比
p1|p2|p3

#---------降维及细胞排序----------
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("05_pseudoTime/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("pseudotime/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("05_pseudoTime/Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("05_pseudoTime/Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("05_pseudoTime/Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("05_pseudoTime/Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plotc <- plot1|plot2|plot3
ggsave("05_pseudoTime/Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("05_pseudoTime/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "05_pseudoTime/pseudotime.csv")

save(mycds, file="./2.results/05_pseudoTime_mycds.rdata")
