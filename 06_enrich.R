library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())
dir.create("06_enrich")

load(file="./2.results/04_subCl_scRNA.rdata") 
load(file="./2.results/05_pseudoTime_mycds.rdata") 

# 基因差异表达分析
#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(scRNA,ident.1 = 0,ident.2 = 1)
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较B_cell和T_cells的差异表达基因
dge.celltype <- FindMarkers(scRNA, ident.1 = 'B_cell', ident.2 = 'T_cells', group.by = 'celltype')
sig_dge.celltype <- subset(dge.celltype, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较拟时State1和State3的差异表达基因
p_data <- subset(pData(mycds),select='State')
scRNAsub <- subset(scRNA, cells=row.names(p_data))
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
dge.State <- FindMarkers(scRNAsub, ident.1 = 1, ident.2 = 3, group.by = 'State')
sig_dge.State <- subset(dge.State, p_val_adj<0.01&abs(avg_log2FC)>1)

#差异基因GO富集分析
ego_ALL <- enrichGO(gene = row.names(sig_dge.celltype),
                    #universe     = row.names(dge.celltype),
                    OrgDb = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'06_enrich/enrichGO.csv')           
ego_CC <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF
ggsave('06_enrich/enrichGO.png', plotc, width = 12,height = 10)


# 差异基kegg富集分析
genelist <- bitr(row.names(sig_dge.celltype), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave("06_enrich/enrichKEGG.png", plot = plotc, width = 12, height = 10)