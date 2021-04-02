suppressMessages(library(Seurat))
suppressMessages(library(MAST))
suppressMessages(library(tidyverse))
suppressMessages(library(pacman))
suppressMessages(library(data.table))
rm(list=ls())
dir.create("X01_QC")
dir.create('./download/01_base_file/', recursive=T)
dir.create('processData')
##----------创建seurat对象------------
inPath = "./download/01_base_file/filtered_feature_bc_matrix"
scRNA.counts <- Read10X(data.dir = inPath)    
try({scRNA = CreateSeuratObject(scRNA.counts[['Gene Expression']])},silent=T)
if(!exists('scRNA')){
  scRNA = CreateSeuratObject(scRNA.counts,project = "SeuratProject")
}
# table(scRNA@meta.data$orig.ident) #查看样本的细胞数量

##--------------计算质控指标-------------
#计算细胞中线粒体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)

col.num <- length(levels(scRNA@active.ident))
ind = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB")
violin <- VlnPlot(scRNA,
                  features = ind, 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) 
ggsave("X01_QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("X01_QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("X01_QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("X01_QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)


#----------设置质控标准-------------
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", 
        "maxGene=4000", "pctMT=20"))
minGene=500
maxGene=4000
pctMT=15

#-------------数据质控------------
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
violin <-VlnPlot(scRNA,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) 
ggsave("X01_QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("X01_QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

#----------------保存中间结果-----------------
save(scRNA, file="./processData/01_QC_scRNA.rdata")
