library(XICOR)
library(Seurat)
library(ggplot2)
library(pbapply)
library(patchwork)
library(reshape2)
library(dplyr)
library(igraph)
library(ggraph)
library(cowplot)
source('code/XICOR_mod.R')
source('code/scripts.R')

#---------------------
# Same 500-gene set used in RealDataCutoffSelection_ver2.R (identical seed/source)
#---------------------
set.seed(0)
commonGenes <- sample(read.csv('data/TcellCommonGenes.csv', row.names = 1)$x, 1000)
commonGenes <- commonGenes[!c(startsWith(commonGenes, 'RPL') | startsWith(commonGenes, 'RPS') | startsWith(commonGenes, 'MT-') | startsWith(commonGenes, 'ENSG'))]
nCells <- 1500

#---------------------
# Load TFLink database (directed ground truth)
#---------------------
TFLink <- read.csv('data/TFTargetDBs/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz', sep = '\t')
TFLinkTFTarget <- TFLink[c('Name.TF', 'Name.Target')] %>%
  mutate('TF.Target' = paste0(Name.TF, "_", Name.Target))
TFLinkTFTarget$database <- 'TFLink'

#---------------------
# Load STRING database (undirected ground truth)
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
# Cutoff for dependence, learned from GSE239592 (transferred, not re-estimated
# on GSM8331607) -- matches the cutoff-selection script's Dataset 1 procedure
#---------------------
soPBMC <- readRDS('data/PBMC/PBMC_Tcells_merged.rds')
so1 <- subset(soPBMC, subset = orig.dataset == 'GSE239592')
soFilt1 <- so1[commonGenes, ]
soFilt1

dfCount1 <- as.data.frame(t(as.matrix(GetAssayData(soFilt1, layer = 'data'))))
set.seed(0)

STRINGedglist_sub1 <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt1) &
                                      STRINGedglist$protein2 %in% rownames(soFilt1) &
                                      STRINGedglist$combined_score > 950, ]
STRINGadj1 <- matrix(0, nrow(soFilt1), nrow(soFilt1))
colnames(STRINGadj1) <- rownames(soFilt1)
rownames(STRINGadj1) <- rownames(soFilt1)
STRINGadj1[as.matrix(STRINGedglist_sub1[, 1:2])] <- 1
STRINGadj1_logical <- apply(STRINGadj1, 2, as.logical)

cutoffs_Xi <- CalcCut(dfCount1, niter = 50, nsamp = nCells, algorithm = 'Xi',
                      alpha = 0.05, val_mat = STRINGadj1_logical, symmetric = TRUE, nCores = 8)
cutoffs_Xi
# cutoffs_Xi <- c('a95' = 0.05516454, 'FDR95' = 0.05563193)

#---------------------
# Load GSM8331607 from the merged, Harmony-corrected T-cell object and
# restrict to the same 500 common genes used for cutoff selection
#---------------------
so3 <- subset(soPBMC, subset = orig.dataset == 'GSM8331607')
soFilt3 <- so3[commonGenes, ]
soFilt3

dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt3, layer = 'data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells), ]
dim(dfCount)

## STRING ground truth restricted to GSM8331607's genes, used to validate edges
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt3) &
                                     STRINGedglist$protein2 %in% rownames(soFilt3) &
                                     STRINGedglist$combined_score > 950, ]
STRINGedglist_sub <- STRINGedglist_sub %>%
  mutate('Gene.Gene' = paste0(protein1, '_', protein2))

#---------------------
# Correlation matrices on ALL nCells GSM8331607 cells (no train/test split),
# averaged over iterations
#---------------------
correlation_matrices_list <- list()
n_iters <- 10
for (i in 1:n_iters) {
  set.seed(i)
  current_cor_matrix <- abs(xicor_mpar(dfCount, nCores = 12))
  correlation_matrices_list[[i]] <- current_cor_matrix
  cat(paste0("Completed repetition ", i, " of ", n_iters, "\n"))
}
summed_cor_matrix <- Reduce('+', correlation_matrices_list)
cor_xi_test <- summed_cor_matrix / n_iters
cor_xi_test[is.na(cor_xi_test)] <- 0

diag(cor_xi_test) <- 0
cor_xi_filtered <- cor_xi_test
cor_xi_filtered[pmax(cor_xi_test, t(cor_xi_test)) < cutoffs_Xi['a95']] <- 0

#---------------------
# Active-regulation (directed) test using the fitted cutoff function (10)
#---------------------
ActRegCut <- 0.455 * (nCells) ** -0.215
cor_xi_diff <- cor_xi_filtered - t(cor_xi_filtered)

# generating edgelist
xi_edge <- melt(cor_xi_diff, value.name = 'diffCor')
xi_edge$Var1 <- as.character(xi_edge$Var1)
xi_edge$Var2 <- as.character(xi_edge$Var2)
xi_edge <- xi_edge %>% arrange(-diffCor) %>% filter(abs(diffCor) > ActRegCut)
xi_edge_graph <- xi_edge %>% arrange(-diffCor) %>% filter(diffCor > ActRegCut) %>%
  mutate("pairs" = paste0(Var1, "_", Var2))
xi_edge_graph <- xi_edge_graph %>%
  mutate('Validation' = case_when(
    pairs %in% STRINGedglist_sub$Gene.Gene & pairs %in% TFLinkTFTarget$TF.Target ~ "Both",
    pairs %in% STRINGedglist_sub$Gene.Gene ~ "STRING",
    pairs %in% TFLinkTFTarget$TF.Target ~ "TF Gene",
    TRUE ~ "None"
  )) %>%
  mutate(Validation = factor(Validation, levels = c('None', 'STRING', 'TF Gene', 'Both'), ordered = TRUE)) %>%
  arrange(desc(Validation))
head(xi_edge_graph, 10)

# Known TF-target links, for directionality validation
xi_edgeSub <- xi_edge %>% mutate("pairs" = paste0(Var1, "_", Var2)) %>%
  filter(pairs %in% TFLinkTFTarget$TF.Target)
xi_edgeSub

xi_edgeSub_graph <- xi_edgeSub %>%
  mutate(Var1 = ifelse(diffCor > 0, xi_edgeSub$Var1, stringr::str_split(pairs, "_", simplify = T)[, 2])) %>%
  mutate(Var2 = ifelse(diffCor > 0, xi_edgeSub$Var2, stringr::str_split(pairs, "_", simplify = T)[, 1]))
xi_edgeSub_graph

signs <- table(sign(xi_edgeSub$diffCor))
paste0("Correct TF identification = ", round(signs['1'] / sum(signs) * 100, 2), '%')

#---------------------
# All significant DEPENDENT links (undirected), with the subset that also
# passes the active-regulation test highlighted in red with an arrow
#---------------------
sigAdj <- pmax(cor_xi_filtered, t(cor_xi_filtered)) > 0
diag(sigAdj) <- FALSE
sigPairs <- which(upper.tri(sigAdj) & sigAdj, arr.ind = TRUE)
genes <- rownames(cor_xi_filtered)

sigEdges <- data.frame(
  Var1 = genes[sigPairs[, 1]],
  Var2 = genes[sigPairs[, 2]],
  pairAB = paste0(genes[sigPairs[, 1]], '_', genes[sigPairs[, 2]]),
  pairBA = paste0(genes[sigPairs[, 2]], '_', genes[sigPairs[, 1]])
)

# Which significant pairs are ALSO directed (present in xi_edge_graph)?
directedPairsSet <- xi_edge_graph$pairs
sigEdges$Directed <- sigEdges$pairAB %in% directedPairsSet | sigEdges$pairBA %in% directedPairsSet

# Which significant pairs are validated in the STRING database?
sigEdges$STRING <- sigEdges$pairAB %in% STRINGedglist_sub$Gene.Gene | sigEdges$pairBA %in% STRINGedglist_sub$Gene.Gene

# Single edge category for legend/colouring (priority: Active regulation > STRING > Dependent only)
sigEdges$EdgeCategory <- factor(
  case_when(
    sigEdges$Directed ~ 'Active regulation',
    sigEdges$STRING   ~ 'STRING-validated',
    TRUE              ~ 'Dependent only'
  ),
  levels = c('Dependent only', 'STRING-validated', 'Active regulation')
)

# For directed pairs, orient Var1->Var2 to match the inferred direction
# (xi_edge_graph stores only the significant, correctly-oriented row)
swap <- sigEdges$Directed & !(sigEdges$pairAB %in% directedPairsSet)
tmp <- sigEdges$Var1[swap]
sigEdges$Var1[swap] <- sigEdges$Var2[swap]
sigEdges$Var2[swap] <- tmp

networkObjAll <- graph_from_edgelist(as.matrix(sigEdges[, c('Var1', 'Var2')]), directed = TRUE)
networkObjAll <- set_edge_attr(networkObjAll, name = 'EdgeType',
                               value = ifelse(sigEdges$Directed, 'Active regulation', 'Dependent only'))
networkObjAll <- set_edge_attr(networkObjAll, name = 'STRING',
                               value = ifelse(sigEdges$STRING, 'STRING-validated', 'Not in STRING'))
networkObjAll <- set_edge_attr(networkObjAll, name = 'EdgeCategory',
                               value = factor(as.character(sigEdges$EdgeCategory),
                                              levels = c('Dependent only', 'STRING-validated', 'Active regulation')))

set.seed(0)
layoutAll <- layout_with_fr(networkObjAll)

p3 <- ggraph(networkObjAll, layout = layoutAll) +
  geom_edge_link(aes(colour = EdgeCategory, filter = EdgeCategory != 'Active regulation'),
                 edge_width = 0.4) +
  geom_edge_link(aes(colour = STRING, filter = STRING == 'STRING-validated'),
                 edge_width = 0.4) +
  geom_node_point(size = 1, colour = alpha('magenta4', 0.8), shape = 16) +
  geom_edge_link(aes(colour = EdgeCategory, filter = EdgeCategory == 'Active regulation'),
                 arrow = arrow(length = unit(0.05, 'in'), type = 'closed'),
                 edge_width = 0.4,
                 end_cap = circle(0.04, 'in'), start_cap = circle(0.02, 'in')) +
  scale_edge_colour_manual(name = 'Edge type',
                           values = c('Dependent only' = alpha('grey75', 0.4),
                                      'STRING-validated' = alpha('dodgerblue3', 0.8),
                                      'Active regulation' = 'firebrick'),
                           drop = FALSE) +
  guides(edge_colour = guide_legend(override.aes = list(edge_width = 1.2))) +
  theme_graph()
  # ggtitle('All Dependent Links')

png('figures/Network_GSM8331607_AllSignificant_DirectedHighlight.png', height = 6, width = 7, res = 300, unit = 'in')
p3
dev.off()

#---------------------
# Full directed GRN, edges colored by validation source
#---------------------
networkObj <- graph_from_edgelist(as.matrix(xi_edge_graph[, 1:2]), directed = T)
networkObj <- set_edge_attr(networkObj, name = "Validation", value = xi_edge_graph$Validation)

set.seed(0)
layout <- layout_with_fr(networkObj)

p1 <- ggraph(networkObj, layout = layout) +
  geom_edge_link(arrow = arrow(length = unit(0.05, 'in'), type = 'closed'),
                 edge_width = 0.5, color = c(alpha('firebrick', 0.5)),
                 end_cap = circle(0.04, 'in'), start_cap = circle(0.02, 'in')) +
  geom_node_point(size = 1.5, colour = alpha('magenta4', 0.7), shape = 16) +
  geom_node_text(aes(label = name), colour = 'black', size = 4, repel = T, max.overlaps = 3) +
  theme_graph()
  # ggtitle('Active Regulations')
p1
png('figures/Network_GSM8331607_AllGene.png', height = 4, width = 4, res = 300, unit = 'in')
p1
dev.off()

#---------------------
# Scatter plot with spline fit for a chosen gene pair, to visualize the
# interaction underlying a directed edge. Defaults to the top-ranked
# active-regulation edge from xi_edge_graph -- edit g1/g2 to inspect others.
#---------------------
g1 <- xi_edge_graph$Var1[1]
g2 <- xi_edge_graph$Var2[1]

# genePairDF <- dfCount[, c(g1, g2)]
genePairDF <- apply(dfCount[, c(g1, g2)], 2, rank, ties.method = 'random')

pSpline <- ggplot(genePairDF, aes(x = .data[[g1]], y = .data[[g2]])) +
  geom_point(size = 1, colour = alpha('black', 0.5), shape = 16) +
  geom_smooth(method = 'loess', #formula = y ~ s(x, bs = 'cs'),
              se = TRUE, colour = 'magenta4', fill = alpha('magenta4', 0.2)) +
  theme_classic() +
  xlab(paste0(g1, ' (Ranked expr.)')) +
  ylab(paste0(g2, ' (Ranked expr.)')) +
  ggtitle(paste0(g1, ' \u2192 ', g2, '\n(\u03be = ', round(cor_xi_test[g1, g2], 3),
                 '  vs  \u03be reverse = ', round(cor_xi_test[g2, g1], 3), ')'))
pSplineFlip <- ggplot(genePairDF, aes(x = .data[[g2]], y = .data[[g1]])) +
  geom_point(size = 1, colour = alpha('black', 0.5), shape = 16) +
  geom_smooth(method = 'loess', #formula = y ~ s(x, bs = 'cs'),
              se = TRUE, colour = 'magenta4', fill = alpha('magenta4', 0.2)) +
  theme_classic() +
  xlab(paste0(g2, ' (Ranked expr.)')) +
  ylab(paste0(g1, ' (Ranked expr.)')) +
  ggtitle(NULL)
pSpline + pSplineFlip

png(paste0('figures/Scatter_GSM8331607_', g1, '_', g2, '.png'), height = 4, width = 7, units = 'in', res = 300)
pSpline + pSplineFlip
dev.off()


#---------------------
# Scatterplot matrix of a few selected genes, to visualize the pairwise
#---------------------

pairsDF <- dfCount[,c('SYNE2','S100A4','LTB','CTSW','IL7R',
                                  'HLA-B','B2M','ATXN7','SIK2')]
pairsDF <- apply(pairsDF, 2, rank, ties.method = 'random')

pairsPlot <- pairs(pairsDF, pch = 16, cex = 0.2, col = alpha('black', 0.5),
      lower.panel = panel.smooth, upper.panel = panel.smooth,
      method = 'loess', col.smooth = 'magenta4', lwd = 1.5)
pairsPlot

png(paste0('figures/Scatter_GSM8331607_All.png'), height = 10, width = 10, units = 'in', res = 300)
pairs(pairsDF, pch = 16, cex = 0.2, col = alpha('black', 0.5),
      lower.panel = panel.smooth, upper.panel = panel.smooth,
      method = 'loess', col.smooth = 'magenta4', lwd = 1.5)
dev.off()

