library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(SingleCellExperiment)

# Helper: run scDblFinder on a Seurat object (raw counts), drop doublets
run_dblfinder <- function(so, seed = 42) {
  set.seed(seed)
  sce <- as.SingleCellExperiment(so, assay = "RNA")   # uses counts layer
  sce <- scDblFinder(sce)                              # add samples= if object has multiple captures
  so$scDblFinder.class <- sce$scDblFinder.class
  so$scDblFinder.score <- sce$scDblFinder.score
  message(sprintf("%s: %d doublets / %d cells (%.1f%%)",
                  unique(so$orig.dataset), sum(so$scDblFinder.class == "doublet"),
                  ncol(so), 100 * mean(so$scDblFinder.class == "doublet")))
  subset(so, subset = scDblFinder.class == "singlet")
}

# Data1 - GSE239592
so1 <- readRDS('data/GSE239592_UPMC_scRNA_integrated_combined_reSCT.rds')
so1_sub <- subset(so1, subset = Type == 'CONTROL')
DefaultAssay(so1_sub) <- 'RNA'
so1_sub@assays$SCT <- NULL
so1_sub <- subset(so1_sub, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so1_sub <- so1_sub[rowSums(GetAssayData(so1_sub, layer = 'counts') > 0) > 10, ]
so1_sub <- UpdateSeuratObject(so1_sub)
so1_sub$orig.dataset <- 'GSE239592'
so1_sub <- run_dblfinder(so1_sub)   # NOTE: see sample-column comment below
so1_sub
rm(so1)
gc()

# Data3 - GSM8331607
mat <- Read10X('data/GSM8331607/')
so3 <- CreateSeuratObject(CreateAssay5Object(mat))
so3[["percent.mt"]] <- PercentageFeatureSet(so3, pattern = "^MT-")
so3 <- subset(so3, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so3 <- so3[rowSums(GetAssayData(so3, layer = 'counts') > 0) > 10, ]
so3$orig.dataset <- 'GSM8331607'
so3 <- run_dblfinder(so3)
so3
gc()

# Data4 - GSM7873657
mat <- Read10X('data/GSM7873657/')
so4 <- CreateSeuratObject(CreateAssay5Object(mat))
so4[["percent.mt"]] <- PercentageFeatureSet(so4, pattern = "^MT-")
so4 <- subset(so4, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so4 <- so4[rowSums(GetAssayData(so4, layer = 'counts') > 0) > 10, ]
so4$orig.dataset <- 'GSM7873657'
so4 <- run_dblfinder(so4)
so4
gc()

# Data5 - GSM7873658
mat <- Read10X('data/GSM7873658/')
so5 <- CreateSeuratObject(CreateAssay5Object(mat))
so5[["percent.mt"]] <- PercentageFeatureSet(so5, pattern = "^MT-")
so5 <- subset(so5, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
so5 <- so5[rowSums(GetAssayData(so5, layer = 'counts') > 0) > 10, ]
so5$orig.dataset <- 'GSM7873658'
so5 <- run_dblfinder(so5)
so5
gc()

#---------------------
# Merge PBMC datasets and batch-correct with Harmony
#---------------------
soPBMC <- merge(so1_sub, y = list(so3, so4, so5),
                add.cell.ids = c('GSE239592', 'GSM8331607', 'GSM7873657', 'GSM7873658'),
                project = 'PBMC_Tcells')
saveRDS(soPBMC, 'data/PBMC/PBMC_Tcells_DoubletRemoved.rds', compress = FALSE)
rm(so1_sub, so3, so4, so5); gc()


#---------------------
# Load Doublet removed data and batch-correct with Harmony
#---------------------
soPBMC <- readRDS('data/PBMC/PBMC_Tcells_DoubletRemoved.rds')
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

soPBMC <- RunUMAP(soPBMC, reduction = 'harmony', umap.method = 'uwot2',
                  dims = 1:20, n.neighbors = 20)

## Cell type annotation - now done ONCE on the merged, Harmony-corrected object
soPBMC <- FindNeighbors(soPBMC, reduction = 'harmony', dims = 1:20, k.param = 20)
soPBMC <- FindClusters(soPBMC, res = 1, cluster.name = 'louvain')
table(soPBMC$louvain, soPBMC$orig.dataset)
Idents(soPBMC) <- 'louvain'
markers <- FindAllMarkers(soPBMC, only.pos = TRUE, min.pct = 0.1,
                          logfc.threshold = 1)

DimPlot(soPBMC, group.by = 'orig.dataset')  # check batch mixing post-Harmony
DimPlot(soPBMC, group.by = 'louvain', label = T)
DotPlot(soPBMC, features = c('CD3D', 'CD3E', 'CD3G', 'MKI67', 'CD4', 'FOXP3',
                             'CTLA4', 'CD8A', 'CD8B', 'SELL','IL7R', 'GZMB','CCR7',
                             'MS4A1', 'CD14', 'FCGR3A', 'CD68', 'ITGAX', 'CD1C',
                             'IL3RA', 'NCR1', 'NKG7','NCAM1', 'PPBP'), group.by = 'louvain',
        cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

soPBMC$celltype <- recode_factor(
  soPBMC$louvain,
  '0'  = 'CD4 T naive',
  '5'  = 'CD4 T naive',
  '7'  = 'CD4 T naive',
  '1'  = 'CD4 T memory',
  '2'  = 'CD4 T memory',          # faint FOXP3/CTLA4, not called Treg
  '10' = 'CD8 T naive',        # CD8A/B high, CCR7+, NKG7-
  '6'  = 'CD8 T cytotoxic',
  '12' = 'CD8 T cytotoxic',
  '17' = 'Doublet',    # weak signal, verify
  '20' = 'Proliferating',         # MKI67 high
  '3'  = 'NK',            # NKG7+, FCGR3A+
  '8'  = 'NK',
  '19' = 'NK',
  '15' = 'NK',         # SELL+, FCGR3A low
  '4'  = 'B',
  '11' = 'B',
  '9'  = 'CD14 Mono',
  '13' = 'CD14 Mono',
  '14' = 'CD16 Mono',             # FCGR3A high, CD14 low
  '18' = 'cDC2',                  # CD1C+
  '21' = 'pDC',                   # IL3RA+
  '16' = 'Platelet'               # PPBP high; CD3D/E likely contamination
)

# Coarse lineage (T vs Other), for T cell extraction
t_labels <- c('CD4 T naive', 'CD4 T memory',
              'CD8 T naive', 'CD8 T cytotoxic')
soPBMC$lineage <- ifelse(soPBMC$celltype %in% t_labels, 'T cells', 'Other')

table(soPBMC$celltype, soPBMC$orig.dataset)

saveRDS(soPBMC, 'data/PBMC/PBMC_annotated.rds')

# NOTE: inspect the DotPlot/UMAP above and set the T-cell cluster IDs for THIS run
soPBMC <- readRDS('data/PBMC/PBMC_annotated.rds')

umapPBMC <- DimPlot(soPBMC, cols = hcl.colors(7, 'RdPu')[c(5, 2)], group.by = 'lineage',
                    pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  ggtitle(NULL)
umapPBMC
umapCelltype <- DimPlot(soPBMC, group.by = 'celltype',
                    pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  ggtitle(NULL)
umapCelltype
umapSample <- DimPlot(soPBMC, group.by = 'orig.dataset',
                      pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  ggtitle(NULL)
umapGene <- FeaturePlot(soPBMC, features = 'CD3D', cols = hcl.colors(7, 'Purples')[c(6, 1)],
                      pt.size = 0.01, alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = 'bottom') +
  ggtitle('CD3D')
umapGene
umapSample + umapPBMC + umapCelltype + umapGene

png(paste0('figures/RealData_PBMC.png'), height = 4, width = 15, units = 'in', res = 300)
umapSample + umapPBMC + umapCelltype +
  patchwork::plot_layout(ncol = 3)
dev.off()

png(paste0('figures/RealData_PBMC_CD3D.png'), height = 5, width = 4.5, units = 'in', res = 300)
umapGene; dev.off()

# marker gene dotplot
# Order top to bottom, matching the marker order below
ct_order <- c('CD4 T naive', 'CD4 T memory', 'CD8 T naive', 'CD8 T cytotoxic',
              'Proliferating', 'NK', 'B', 'CD14 Mono', 'CD16 Mono',
              'cDC2', 'pDC', 'Platelet', 'Doublet')
soPBMC$celltype <- factor(soPBMC$celltype, levels = rev(ct_order))
markers <- c('CD3D', 'CD3E', 'CD3G', 'CD4', 'SELL', 'CCR7', 'IL7R',  # CD4 T
             'FOXP3', 'CTLA4',                                        # CD4 memory/Treg-like
             'CD8A', 'CD8B',                                          # CD8 naive
             'NKG7', 'GZMB',                                          # cytotoxic / NK
             'MKI67',                                                 # proliferating
             'NCR1', 'NCAM1',                                         # NK
             'MS4A1',                                                 # B
             'CD14',                                                  # CD14 Mono
             'FCGR3A', 'CD68',                                        # CD16 Mono
             'ITGAX', 'CD1C',                                         # cDC2
             'IL3RA',                                                 # pDC
             'PPBP')                                                  # Platelet

png(paste0('figures/RealData_PBMC_Dotplot.png'), height = 5, width = 10, units = 'in', res = 300)
DotPlot(soPBMC, features = markers, group.by = 'celltype')+
  scale_color_gradientn(colours = hcl.colors(7, 'RdPu')[c(7, 2)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab(NULL)
dev.off()

soTcells <- subset(soPBMC, subset = lineage == 'T cells')
cnt  <- GetAssayData(soTcells, layer = 'counts')
ds   <- soTcells$orig.dataset[colnames(cnt)]
det  <- sapply(unique(ds), function(d) Matrix::rowMeans(cnt[, ds == d, drop = FALSE] > 0))  # genes x datasets
keep <- rownames(det)[apply(det, 1, min) > 0.3]
message(sprintf('Genes detected in > 30%% of cells in all %d datasets: %d / %d',
                ncol(det), length(keep), nrow(det)))
soTcells <- soTcells[keep, ]
soTcells
soTcells$nCount_RNA  <- Matrix::colSums(GetAssayData(soTcells, layer = 'counts'))
soTcells <- subset(soTcells, subset = nCount_RNA > 10)
soTcells

table(soTcells$orig.dataset)  # confirms GSE/GSM provenance is retained per cell

saveRDS(soTcells, 'data/PBMC/PBMC_Tcells_merged.rds')

#---------------------
# find common genes
#---------------------
soTcells <- readRDS('data/PBMC/PBMC_Tcells_merged.rds')

TcellCommonGenes <- rownames(soTcells)  # already common: single merged object, no per-dataset intersection needed
write.csv(TcellCommonGenes, 'data/TcellCommonGenes.csv')

