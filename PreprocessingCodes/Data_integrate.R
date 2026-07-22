library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)

#---------------------
# Load & QC individual PBMC datasets (no per-dataset annotation yet)
#---------------------
# Data1 - GSE239592
so1 <- readRDS('data/GSE239592_UPMC_scRNA_integrated_combined_reSCT.rds')
so1_sub <- subset(so1, subset = Type == 'CONTROL')
DefaultAssay(so1_sub) <- 'RNA'
so1_sub@assays$SCT <- NULL
so1_sub <- subset(so1_sub, subset = nFeature_RNA > 500 &
                    nCount_RNA > 1000 &
                    percent.mt < 10)
so1_sub <- so1_sub[rowSums(GetAssayData(so1_sub, layer = 'counts') > 0) > 10, ]
so1_sub <- UpdateSeuratObject(so1_sub)  # drop stale reductions/assays before merge
so1_sub$orig.dataset <- 'GSE239592'
so1_sub
rm(so1)

# Data3 - GSM8331607
mat <- Read10X('data/GSM8331607/')
so3 <- CreateSeuratObject(CreateAssay5Object(mat))
so3[["percent.mt"]] <- PercentageFeatureSet(so3, pattern = "^MT-")
so3 <- subset(so3, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so3 <- so3[rowSums(GetAssayData(so3, layer = 'counts') > 0) > 10, ]
so3$orig.dataset <- 'GSM8331607'
so3

# Data4 - GSM7873657
mat <- Read10X('data/GSM7873657/')
so4 <- CreateSeuratObject(CreateAssay5Object(mat))
so4[["percent.mt"]] <- PercentageFeatureSet(so4, pattern = "^MT-")
so4 <- subset(so4, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so4 <- so4[rowSums(GetAssayData(so4, layer = 'counts') > 0) > 10, ]
so4$orig.dataset <- 'GSM7873657'
so4

# Data5 - GSM7873658
mat <- Read10X('data/GSM7873658/')
so5 <- CreateSeuratObject(CreateAssay5Object(mat))
so5[["percent.mt"]] <- PercentageFeatureSet(so5, pattern = "^MT-")
so5 <- subset(so5, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so5 <- so5[rowSums(GetAssayData(so5, layer = 'counts') > 0) > 10, ]
so5$orig.dataset <- 'GSM7873658'
so5

gc()
#---------------------
# Merge PBMC datasets and batch-correct with Harmony
#---------------------
soPBMC <- merge(so1_sub, y = list(so3, so4, so5),
                add.cell.ids = c('GSE239592', 'GSM8331607', 'GSM7873657', 'GSM7873658'),
                project = 'PBMC_Tcells')

# rm(so1_sub, so3, so4, so5); gc()
soPBMC
soPBMC <- NormalizeData(soPBMC, normalization.method = "LogNormalize", scale.factor = 10000)
soPBMC <- FindVariableFeatures(soPBMC)
soPBMC <- ScaleData(soPBMC)
soPBMC <- RunPCA(soPBMC, npcs = 30)

# Harmony integration across the 4 samples (Seurat v5 layer-aware API)
soPBMC <- IntegrateLayers(object = soPBMC, method = HarmonyIntegration,
                          orig.reduction = 'pca', new.reduction = 'harmony',
                          max.iter.harmony = 20L)
soPBMC <- JoinLayers(soPBMC)  # collapse per-sample count/data layers now that integration is done

soPBMC <- RunUMAP(soPBMC, reduction = 'harmony', dims = 1:20, n.neighbors = 20)

## Cell type annotation - now done ONCE on the merged, Harmony-corrected object
soPBMC <- FindNeighbors(soPBMC, reduction = 'harmony', dims = 1:20, k.param = 20)
soPBMC <- FindClusters(soPBMC, res = 1, cluster.name = 'louvain')
DimPlot(soPBMC, group.by = 'orig.dataset')  # check batch mixing post-Harmony
DimPlot(soPBMC, group.by = 'louvain', label = T)
DotPlot(soPBMC, features = c('CD3D', 'CD3E', 'CD3G','CD4', 'CD8A', 'CD8B',
                             'MS4A1', 'CD14', 'ITGAX','NCR1'), group.by = 'louvain')

# NOTE: inspect the DotPlot/UMAP above and set the T-cell cluster IDs for THIS run
Tcell_clusters <- c('0', '1', '2', '5', '6', '8', '9', '12')  # <- update after inspection
ct <- ifelse(soPBMC$louvain %in% Tcell_clusters, 'T cells', 'Other')
soPBMC$celltype <- factor(ct)
table(soPBMC$celltype, soPBMC$orig.dataset)

umapPBMC <- DimPlot(soPBMC, cols = hcl.colors(7, 'RdPu')[c(5, 2)], group.by = 'celltype',
                    pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  ggtitle('Harmony')
umapSample <- DimPlot(soPBMC, group.by = 'orig.dataset',
                      pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  ggtitle('Dataset')
umapGene <- FeaturePlot(soPBMC, features = 'CD3D', cols = hcl.colors(7, 'RdPu')[c(5, 1)],
                      pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  ggtitle('CD3D')
umapSample + umapPBMC + umapGene

png(paste0('figures/RealData_PBMC.png'), height = 4, width = 14, units = 'in', res = 300)
umapSample + umapPBMC + umapGene + patchwork::plot_layout(ncol = 3)
dev.off()

soTcells <- subset(soPBMC, subset = celltype == 'T cells')
soTcells <- soTcells[rowSums(GetAssayData(soTcells, layer = 'counts') > 0) > ncol(soTcells) * 0.25, ]
soTcells
table(soTcells$orig.dataset)  # confirms GSE/GSM provenance is retained per cell

saveRDS(soTcells, 'data/PBMC/PBMC_Tcells_merged.rds')

#---------------------
# find common genes
#---------------------
soTcells <- readRDS('data/PBMC/PBMC_Tcells_merged.rds')

TcellCommonGenes <- rownames(soTcells)  # already common: single merged object, no per-dataset intersection needed
write.csv(TcellCommonGenes, 'data/TcellCommonGenes.csv')
