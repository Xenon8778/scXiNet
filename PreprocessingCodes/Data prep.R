library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
source('code/scTypeAnnotation_v2.r')

#---------------------
# Load all data sets
#---------------------
# Data1
so1 <- readRDS('data/GSE239592_UPMC_scRNA_integrated_combined_reSCT.rds')
DefaultAssay(so1) <- 'RNA'
so1@assays$SCT <- NULL

so1 <- subset(so1, subset=Type=='CONTROL')
so1
DimPlot(so1, group.by = 'cellType.SCT')

so1 <- NormalizeData(so1, normalization.method = "LogNormalize",
                     scale.factor = 10000)
so1 <- FindVariableFeatures(so1)
so1 <- ScaleData(so1)
so1 <- RunPCA(so1)
so1 <- RunUMAP(so1, dims = 1:15)

## Cell type annotation
so1 <- FindNeighbors(so1)
so1 <- FindClusters(so1, res=0.8, cluster.name = 'louvain')
DimPlot(so1, group.by = 'louvain')
DotPlot(so1, features = c('CD3D','CD3E','CD3G'), group.by = 'louvain')

ct <- so1$louvain
ct <- recode_factor(ct, '0'='T cells', '1'='T cells', '2'='T cells',
                    '5'='T cells', '7'='T cells', '9' = 'T cells',
                    '13' = 'T cells', '14' = 'T cells')
ct <- as.character(ct)
ct <- sapply(ct, function(x) ifelse(x == 'T cells', 'T cells', 'Other'),
             USE.NAMES = F)
so1$celltype <- factor(ct)

umap1 <- DimPlot(so1, cols = hcl.colors(7, 'RdPu')[c(5,2)], group.by = 'celltype')+NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  ggtitle('GSE239592')
umap1

soTcells <- subset(so1, subset = celltype == 'T cells')
soTcells <- soTcells[rowSums(GetAssayData(soTcells, layer='count')>0) > ncol(soTcells)*0.20,]
soTcells

saveRDS(soTcells, 'data/PBMC/GSE239592_Tcells.rds')

# Data3
mat <- Read10X('data/GSM8331607/')
so3 <- CreateSeuratObject(CreateAssayObject(mat))
so3[["percent.mt"]] <- PercentageFeatureSet(so3, pattern = "^MT-")
so3 <- subset(so3, subset = nCount_RNA > 1000 & percent.mt < 10)
so3 <- NormalizeData(so3, normalization.method = "LogNormalize",
                     scale.factor = 10000)
so3 <- FindVariableFeatures(so3)
so3 <- ScaleData(so3)
so3 <- RunPCA(so3)
so3 <- RunUMAP(so3, dims = 1:15)

## Cell type annotation
so3 <- FindNeighbors(so3)
so3 <- FindClusters(so3, res=0.8, cluster.name = 'louvain')
DimPlot(so3, group.by = 'louvain')
DotPlot(so3, features = c('CD3D','CD3E','CD3G'), group.by = 'louvain')
ct <- so3$louvain
ct <- recode_factor(ct, '2'='T cells', '4'='T cells',
                    '5'='T cells', '6'='T cells', '9' = 'T cells',
                    '13' = 'T cells')
ct <- as.character(ct)
ct <- sapply(ct, function(x) ifelse(x == 'T cells', 'T cells', 'Other'),
             USE.NAMES = F)
so3$celltype <- factor(ct)
table(so3$celltype)

umap3 <- DimPlot(so3, cols = hcl.colors(7, 'RdPu')[c(5,2)], group.by = 'celltype')+NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  ggtitle('GSM8331607')
umap3

soTcells <- subset(so3, subset = celltype == 'T cells')
soTcells <- soTcells[rowSums(GetAssayData(soTcells, layer='count')>0) > ncol(soTcells)*0.20,]
soTcells

saveRDS(soTcells, 'data/PBMC/GSM8331607_Tcells.rds')

# Data4
mat <- Read10X('data/GSM7873657/')
so4 <- CreateSeuratObject(CreateAssayObject(mat))
so4[["percent.mt"]] <- PercentageFeatureSet(so4, pattern = "^MT-")
so4 <- subset(so4, subset = nCount_RNA > 1000 & percent.mt < 10)
so4 <- NormalizeData(so4, normalization.method = "LogNormalize",
                     scale.factor = 10000)
so4 <- FindVariableFeatures(so4)
so4 <- ScaleData(so4)
so4 <- RunPCA(so4)
so4 <- RunUMAP(so4, dims = 1:15)

## Cell type annotation
so4 <- FindNeighbors(so4)
so4 <- FindClusters(so4, res=0.8, cluster.name = 'louvain')
DimPlot(so4, group.by = 'louvain')
DotPlot(so4, features = c('CD3D','CD3E','CD3G'), group.by = 'louvain')
ct <- so4$louvain
ct <- recode_factor(ct, '0'='T cells', '2'='T cells', '4'='T cells',
                    '5'='T cells', '6'='T cells', '7' = 'T cells',
                    '8' = 'T cells','9' = 'T cells','10'='T cells',
                    '14' = 'T cells')
ct <- as.character(ct)
ct <- sapply(ct, function(x) ifelse(x == 'T cells', 'T cells', 'Other'),
             USE.NAMES = F)
so4$celltype <- factor(ct)

umap4 <- DimPlot(so4, cols = hcl.colors(7, 'RdPu')[c(5,2)], group.by = 'celltype')+NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  ggtitle('GSM7873657')
umap4

soTcells <- subset(so4, subset = celltype == 'T cells')
soTcells <- soTcells[rowSums(GetAssayData(soTcells, layer='count')>0) > ncol(soTcells)*0.20,]
soTcells

saveRDS(soTcells, 'data/PBMC/GSM7873657_Tcells.rds')

# Data5
mat <- Read10X('data/GSM7873658/')
so5 <- CreateSeuratObject(CreateAssayObject(mat))
so5[["percent.mt"]] <- PercentageFeatureSet(so5, pattern = "^MT-")
so5 <- subset(so5, subset = nCount_RNA > 1000 & percent.mt < 10)
so5 <- NormalizeData(so5, normalization.method = "LogNormalize",
                     scale.factor = 10000)
so5 <- FindVariableFeatures(so5)
so5 <- ScaleData(so5)
so5 <- RunPCA(so5)
so5 <- RunUMAP(so5, dims = 1:15)

## Cell type annotation
so5 <- FindNeighbors(so5)
so5 <- FindClusters(so5, res=0.8, cluster.name = 'louvain')
DimPlot(so5, group.by = 'louvain')
DotPlot(so5, features = c('CD3D','CD3E','CD3G'), group.by = 'louvain')

ct <- so5$louvain
ct <- recode_factor(ct, '1'='T cells', '2'='T cells', '3'='T cells',
                    '4'='T cells', '5'='T cells', '8' = 'T cells',
                    '11'='T cells')
ct <- as.character(ct)
ct <- sapply(ct, function(x) ifelse(x == 'T cells', 'T cells', 'Other'),
             USE.NAMES = F)
so5$celltype <- factor(ct)

umap5 <- DimPlot(so5, cols = hcl.colors(7, 'RdPu')[c(5,2)], group.by = 'celltype')+NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  ggtitle('GSM7873658')
umap5

soTcells <- subset(so5, subset = celltype == 'T cells')
soTcells <- soTcells[rowSums(GetAssayData(soTcells, layer='count')>0) > ncol(soTcells)*0.20,]
soTcells

saveRDS(soTcells, 'data/PBMC/GSM7873658_Tcells.rds')

# Data8
so8 <- readRDS('data/GSE284419_small_intestine/GSM8683884_seuratObject.CHGA-Venus.duodenum.RDS')
DefaultAssay(so8) <- 'RNA'
so8@assays$SCT <- NULL
so8[["percent.mt"]] <- PercentageFeatureSet(so8, pattern = "^MT-")
so8 <- subset(so8, subset = nCount_RNA > 1000 & percent.mt < 10)
so8 <- NormalizeData(so8, normalization.method = "LogNormalize",
                     scale.factor = 10000)
so8 <- so8[rowSums(GetAssayData(so8, layer='count')>0) > ncol(so8)*0.20,]
so8
so8 <- FindVariableFeatures(so8)
so8 <- ScaleData(so8)
so8 <- RunPCA(so8)
so8 <- RunUMAP(so8, dims = 1:10)
Idents(so8) <- '0'
umap8 <- DimPlot(so8, cols = hcl.colors(7, 'RdPu')[2], group.by = NULL)+NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  ggtitle('GSM8683884')
saveRDS(so8, 'data/Intestine/GSM8683884_Intestine.rds')

# Data9
so9 <- readRDS('data/GSE284419_small_intestine/GSM8683885_seuratObject.CHGA-Venus.ileum.RDS')
DefaultAssay(so9) <- 'RNA'
so9@assays$SCT <- NULL
so9[["percent.mt"]] <- PercentageFeatureSet(so9, pattern = "^MT-")
so9 <- subset(so9, subset = nCount_RNA > 1000 & percent.mt < 10)
so9 <- NormalizeData(so9, normalization.method = "LogNormalize",
                     scale.factor = 10000)
so9 <- so9[rowSums(GetAssayData(so9, layer='count')>0) > ncol(so9)*0.20,]
so9
so9 <- FindVariableFeatures(so9)
so9 <- ScaleData(so9)
so9 <- RunPCA(so9)
so9 <- RunUMAP(so9, dims = 1:10)
Idents(so9) <- '0'
umap9 <- DimPlot(so9, cols = hcl.colors(7, 'RdPu')[2], group.by = NULL)+NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  ggtitle('GSM8683885')
saveRDS(so9, 'data/Intestine/GSM8683885_Intestine.rds')

message('Loaded samples removed...')
rm(so1,so3,so4,so5,so8,so9)

#---------------------
# find 300 common genes
#---------------------
so1 <- readRDS('data/PBMC/GSE239592_Tcells.rds')
so3 <- readRDS('data/PBMC/GSM8331607_Tcells.rds')
so4 <- readRDS('data/PBMC/GSM7873657_Tcells.rds')
so5 <- readRDS('data/PBMC/GSM7873658_Tcells.rds')
so8 <- readRDS('data/Intestine/GSM8683884_Intestine.rds')
so9 <- readRDS('data/Intestine/GSM8683885_Intestine.rds')

TcellCommonGenes <- Reduce(intersect, list(rownames(so1),rownames(so6),
                                         rownames(so4),rownames(so5)))
write.csv(TcellCommonGenes,'data/TcellCommonGenes.csv')

intersectGenes <- Reduce(intersect, list(rownames(so1),rownames(so6),
                                         rownames(so4),rownames(so5),
                                         rownames(so8),rownames(so9)))

# Selecting random 300 common genes
set.seed(0)
commonGenes <- intersectGenes[sample(1:length(intersectGenes), 300)]
write.csv(commonGenes,'data/300CommonGenes.csv')

png(paste0('figures/RealData_PBMC.png'),height = 3, width = 12, units = 'in', res = 300)
umap1+umap3+umap4+umap5+patchwork::plot_layout(ncol = 4)
dev.off()

png(paste0('figures/RealData_Intestine.png'),height = 3, width = 6, units = 'in', res = 300)
umap8+umap9+plot_layout(ncol = 2)
dev.off()
message('Plots made...')

message('Loaded samples removed...')
rm(so1,so3,so4,so5,so8,so9)


