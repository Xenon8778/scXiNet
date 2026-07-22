library(XICOR)
library(Seurat)
library(ggplot2)
library(pbapply)
library(patchwork)
source('code/XICOR_mod.R')
source('code/scripts.R')

set.seed(0)
commonGenes <- sample(read.csv('data/TcellCommonGenes.csv', row.names = 1)$x, 1000)
write.csv(commonGenes, 'data/TcellCommonGenes_1000.csv', row.names = F)
nCells <- 1500

#---------------------
# Load STRING database
#---------------------
STRINGedglist <- read.table('data/STRING/9606.protein.links.v12.0.txt.gz',
                            sep = ' ', header = T)
ProteinInfo <- read.csv('data/STRING/9606.protein.info.v12.0.txt.gz',
                        sep = '\t', header = T)
ProteinMap <- c(ProteinInfo$preferred_name)
names(ProteinMap) <- ProteinInfo$X.string_protein_id

STRINGedglist$protein1 <- ProteinMap[STRINGedglist$protein1]
STRINGedglist$protein2 <- ProteinMap[STRINGedglist$protein2]

#---------------------
# Helper: build STRING adjacency for a given filtered object
#---------------------
buildSTRINGadj <- function(soFilt) {
  sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt) &
                         STRINGedglist$protein2 %in% rownames(soFilt) &
                         STRINGedglist$combined_score > 950, ]
  adj <- matrix(0, nrow(soFilt), nrow(soFilt))
  colnames(adj) <- rownames(soFilt)
  rownames(adj) <- rownames(soFilt)
  adj[as.matrix(sub[, 1:2])] <- 1
  adj_logical <- apply(adj, 2, as.logical)
  list(adj = adj, adj_logical = adj_logical)
}

#---------------------
# Helper: sample nCells cells (no train/test split)
#---------------------
sampleCells <- function(soFilt, nCells, seed = 0) {
  dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer = 'data'))))
  set.seed(seed)
  dfCount[sample(rownames(dfCount), nCells), ]
}

#---------------------
# Helper: FPR / F1 curve vs cutoff for a dataset, given a reference cutoff to mark
#---------------------
getPrecision <- function(corData, c, STRINGadj) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  c(cmI$byClass, cmI$overall)
}

evalCurve <- function(cor_xi, STRINGadj, cutoff_seq = seq(0.01, 0.99, 0.01)) {
  Precisions <- as.data.frame(t(pbsapply(cutoff_seq, FUN = getPrecision,
                                         corData = cor_xi, STRINGadj = STRINGadj,
                                         USE.NAMES = T)))
  Precisions$cutoff <- cutoff_seq
  Precisions
}

#---------------------
# Load merged, Harmony-corrected T-cell object (all 4 PBMC samples)
# NOTE: Harmony only corrects the PCA/UMAP embedding used for clustering;
# it does not alter the log-normalized 'data' layer used below for Xi
# correlation, so per-sample Xi values are unaffected by the batch correction.
# Samples are still distinguished via the retained 'orig.dataset' metadata.
#---------------------
soPBMC <- readRDS('data/PBMC/PBMC_Tcells_merged.rds')

#---------------------
# Dataset 1 (GSE239592): learn the cutoff using ALL 2000 cells
# No train/test split -- CalcCut already bootstraps internally (Algorithm 1)
#---------------------
so1 <- subset(soPBMC, subset = orig.dataset == 'GSE239592')
soFilt1 <- so1[commonGenes, ]

s1 <- buildSTRINGadj(soFilt1)
sum(s1[["adj"]])
dfCount1 <- sampleCells(soFilt1, nCells)

cutoffs_Xi <- CalcCut(dfCount1, niter = 50, nsamp = nCells, algorithm = 'Xi',
                      alpha = 0.05, val_mat = s1$adj_logical, symmetric = TRUE, nCores = 12)
cutoffs_Xi

# Correlation matrix on the full 2000 cells, for visualizing dataset 1's own FPR/F1 curve
cor_xi_1 <- abs(xicor_mpar(dfCount1, nCores = 8))
cor_xi_1[is.na(cor_xi_1)] <- 0
cor_xi_1 <- pmax(cor_xi_1, t(cor_xi_1))
diag(cor_xi_1) <- 0

Precisions1 <- evalCurve(cor_xi_1, s1$adj)

p1 <- ggplot(Precisions1, aes(x = cutoff)) +
  geom_line(aes(y = 1 - Specificity), linewidth = 1, color = 'black') +
  geom_vline(xintercept = cutoffs_Xi['a95'], color = 'firebrick', linetype = 'dashed', linewidth = 0.7) +
  geom_hline(yintercept = 0.05, color = 'firebrick', linewidth = 0.7) +
  theme_classic() + xlab('Cutoff') + ylab('FPR') +
  ggtitle('GSE239592')

q1 <- ggplot(Precisions1, aes(x = cutoff)) +
  geom_line(aes(y = F1), linewidth = 1, color = 'black') +
  geom_vline(xintercept = cutoffs_Xi['a95'], color = 'firebrick', linewidth = 0.7) +
  theme_classic() + xlab('Cutoff') +
  ggtitle('GSE239592')

#---------------------
# Datasets 3, 4, 5: pure evaluation -- apply dataset-1 cutoff to ALL 2000 cells
# (no cutoff learning, no split -- these are held-out test datasets)
# Subset by orig.dataset from the merged object rather than reading separate files
#---------------------
targetDatasets <- c('GSM8331607', 'GSM7873657', 'GSM7873658')

evalTarget <- function(dataset_id, soPBMC, nCells, cutoffs_Xi, title) {
  so <- subset(soPBMC, subset = orig.dataset == dataset_id)
  soFilt <- so[commonGenes, ]
  s <- buildSTRINGadj(soFilt)
  dfCount <- sampleCells(soFilt, nCells)

  cor_xi <- abs(xicor_mpar(dfCount, nCores = 8))
  cor_xi[is.na(cor_xi)] <- 0
  cor_xi <- pmax(cor_xi, t(cor_xi))
  diag(cor_xi) <- 0

  cor_xi_filtered <- apply(cor_xi > cutoffs_Xi['a95'], c(1, 2), as.integer)
  cm_xi <- confusionMatrix(table(cor_xi_filtered, s$adj), positive = "1")
  print(cm_xi)

  Precisions <- evalCurve(cor_xi, s$adj)

  p <- ggplot(Precisions, aes(x = cutoff)) +
    geom_line(aes(y = 1 - Specificity), linewidth = 1, color = 'black') +
    geom_vline(xintercept = cutoffs_Xi['a95'], color = 'firebrick', linetype = 'dashed', linewidth = 0.7) +
    geom_hline(yintercept = 0.05, color = 'firebrick', linewidth = 0.7) +
    theme_classic() + xlab('Cutoff') + ylab('FPR') +
    ggtitle(title)

  q <- ggplot(Precisions, aes(x = cutoff)) +
    geom_line(aes(y = F1), linewidth = 1, color = 'black') +
    geom_vline(xintercept = cutoffs_Xi['a95'], color = 'firebrick', linewidth = 0.7) +
    theme_classic() + xlab('Cutoff') +
    ggtitle(title)

  list(p = p, q = q, cm = cm_xi, precisions = Precisions)
}

res2 <- evalTarget('GSM8331607', soPBMC, nCells, cutoffs_Xi, 'GSM8331607')
res3 <- evalTarget('GSM7873657', soPBMC, nCells, cutoffs_Xi, 'GSM7873657')
res4 <- evalTarget('GSM7873658', soPBMC, nCells, cutoffs_Xi, 'GSM7873658')

p2 <- res2$p; q2 <- res2$q
p3 <- res3$p; q3 <- res3$q
p4 <- res4$p; q4 <- res4$q
p1 + p2 + p3 + p4 + plot_layout(ncol = 4) &
  ylim(0, 0.4) & xlim(0, 0.3)

#---------------------
# Comparison plots
#---------------------
png('figures/RealData_FPR_v2.png', height = 3, width = 10, units = 'in', res = 300)
p1 + p2 + p3 + p4 + plot_layout(ncol = 4) &
  ylim(0, 0.4) & xlim(0, 0.3)
dev.off()

png('figures/RealData_F1Score_v2.png', height = 3, width = 10, units = 'in', res = 300)
q1 + q2 + q3 + q4 + plot_layout(ncol = 4) &
  xlim(0, 0.8)
dev.off()


#---------------------
# Barplot of the observed FPR (at cutoffs_Xi['a95']) across the three held-out datasets
#---------------------
fprDF <- data.frame(
  dataset = c('GSM8331607', 'GSM7873657', 'GSM7873658'),
  FPR = c(1 - res2$cm$byClass['Specificity'],
          1 - res3$cm$byClass['Specificity'],
          1 - res4$cm$byClass['Specificity'])
)

fprBox <- ggplot(fprDF, aes(x = 'Held-out datasets', y = FPR)) +
  geom_col(aes(fill = dataset), color='black', position = position_dodge(1), width = 0.8,) +
  geom_hline(yintercept = 0.05, color = 'firebrick', linetype = 'dashed', linewidth = 0.7) +
  theme_classic() + labs(x = NULL , y='FPR', fill = 'Dataset') +
  lims(y = c(0, 0.25)) +
  scale_fill_manual(values = scales::hue_pal()(4)[2:4]) +
  ggtitle('FPR\n(held-out datasets)')
fprBox

png('figures/RealData_FPR_boxplot.png', height = 3, width = 3, units = 'in', res = 300)
fprBox
dev.off()
