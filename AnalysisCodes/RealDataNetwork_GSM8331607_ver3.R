library(XICOR)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(patchwork)
library(ggh4x)
source('code/realDataSetup.R')   # gene panel, STRING GT, sampleCells, computeXi, alpha, nCells

targetDataset <- 'GSM8331607'

# ---- Test choices ----
dep_crit   <- 'FPR'   # dependency cutoff from GSE239592: 'FPR' (FPR) or 'FDR' (FDR95)
asy_crit   <- 'FDR'   # asymmetric-association cutoff function: 'FPR' (FPR) or 'FDR'.
# real dependent pairs are likely mostly monotonic -> use FPR by default.

# ---- Asymmetric-association cutoff ----
asy_cut_method <- 'learned'   # 'learned' = Algorithm 2 on the TARGET dataset at n = nCells,
#             N0 = dependent CORUM co-complex pairs; FPR criterion, asy_crit is ignored
# 'formula' = simulation-fitted c1 * n^-c2 (uses asy_crit)
learn_cores    <- 8          # workers for LearnCutoff
learn_B        <- 50          # Algorithm 2 resampling iterations

# ---- Statistic: ChatterjeeKNN (revised xi, Lin & Han 2023). MUST match RealDataCutoffSelection_ver3.R ----
xi_alg <- 'Xi_KNN'
xi_M   <- 5                          # right nearest neighbours (same as cutoff-selection script)
tag    <- paste0(xi_alg, '_ver3')     # -> Xi_KNN_ver3, same tag as the cutoff-selection outputs

# ---- Reproduce original script (defaults) vs reuse cutoff-selection outputs ----
use_saved_cutoff <- TRUE  # read dependency cutoff from RealDataCutoffSelection_ver3.R (same scheme)
ref_all_cells    <- TRUE  # fallback CalcCut: pool = ALL GSE239592 cells, subsample nCells w/o replacement

#---------------------
# Ground truth
#   STRING three-class GT (shared with cutoff selection) -> undirected validation
#   TFLink TF -> target                                   -> direction validation
#---------------------
gt <- LoadGT_STRING(commonGenes)

TFLink <- read.csv('data/TFTargetDBs/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz',
                   sep = '\t')
tf_pairs <- unique(paste0(TFLink$Name.TF, '_', TFLink$Name.Target))

# N0 for Algorithm 2 (no active regulation): subunits of the same CORUM complex
# (co-regulated, no directed regulation), with no TFLink link in either direction.
corum <- read.delim('data/corum_humanComplexes.txt')        # check column name below in your file
cx    <- lapply(strsplit(corum$subunits_gene_name, ';'), function(s) intersect(unique(trimws(s)), commonGenes))
cx    <- cx[lengths(cx) >= 2]
cx_pairs <- unique(unlist(lapply(cx, function(s) {
  p <- combn(sort(s), 2); paste0(p[1, ], '_', p[2, ])
})))
message(sprintf('CORUM: %d complexes on panel, %d co-complex pairs', length(cx), length(cx_pairs)))

asyNull <- function(g1, g2) {
  a <- pmin(g1, g2); b <- pmax(g1, g2)                       # same sorted order as cx_pairs
  paste0(a, '_', b) %in% cx_pairs &
    !(paste0(g1, '_', g2) %in% tf_pairs | paste0(g2, '_', g1) %in% tf_pairs)
}

#---------------------
# Dependency cutoff: transferred from GSE239592 (RealDataCutoffSelection_ver3.R).
# Re-estimated here only if that output is missing.
#---------------------
cut_file <- paste0('Results/RealData_Cutoffs_GSE239592_', tag, '.csv')
soPBMC <- readRDS('data/PBMC/PBMC_Tcells_merged.rds')
table(soPBMC$orig.dataset)
getCounts <- function(dataset_id, n = nCells) {
  so <- subset(soPBMC, subset = orig.dataset == dataset_id)
  sampleCells(so[commonGenes, ], n)
}

#---------------------
# Dataset diagnostics: how the typical gene on the panel looks in each dataset.
# Raw counts are the scale of the sensitivity simulations (mean_expression, dispersion, dropout);
# the log-normalised 'data' layer is what xi is computed on (sampleCells).
#   zero_frac : fraction of cells with 0 (observed dropout; includes sampling zeros, unlike the
#               simulation's dropout_rate, which is added on top of Poisson/NB zeros)
#   nb_disp   : per-gene method-of-moments NB overdispersion, (var - mean) / mean^2 (raw scale only;
#               crude: no size-factor correction, so inflated by library-size variation)
#   cv2_gene_means : CV^2 of gene means across genes = analogue of the simulation's 'dispersion'
#                    (gene means ~ Gamma(shape = 1/dispersion) -> CV^2 = dispersion)
# All cells of each dataset are used (not the nCells subsample).
#---------------------
diag_datasets <- c(refDataset, 'GSM7873658', 'GSM7873657', 'GSM8331607')
geneStats <- function(m) {                                   # m = genes x cells (sparse ok)
  mu <- Matrix::rowMeans(m)
  v  <- Matrix::rowMeans(m^2) - mu^2
  data.frame(gene = rownames(m), zero_frac = 1 - Matrix::rowSums(m > 0) / ncol(m),
             mean = mu, sd = sqrt(pmax(v, 0)),
             nb_disp = ifelse(mu > 0, pmax(v - mu, 0) / mu^2, NA_real_))
}
diag_genes <- bind_rows(lapply(diag_datasets, function(ds) {
  so <- subset(soPBMC, subset = orig.dataset == ds)[commonGenes, ]
  so <- tryCatch(SeuratObject::JoinLayers(so), error = function(e) so)   # Seurat v5 split layers
  bind_rows(
    data.frame(Dataset = ds, Scale = 'Raw counts', nCellsTotal = ncol(so),
               geneStats(GetAssayData(so, layer = 'counts'))),
    data.frame(Dataset = ds, Scale = 'Log-normalised', nCellsTotal = ncol(so),
               geneStats(GetAssayData(so, layer = 'data'))))
}))
diag_genes$Dataset <- factor(diag_genes$Dataset, levels = diag_datasets)
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
diag_summary <- diag_genes %>%
  group_by(Dataset, Scale, nCellsTotal) %>%
  summarise(nGenes = n(),
            zero_frac_median = median(zero_frac), zero_frac_q25 = q(zero_frac, .25), zero_frac_q75 = q(zero_frac, .75),
            mean_median = median(mean), mean_q25 = q(mean, .25), mean_q75 = q(mean, .75),
            sd_median = median(sd), sd_q25 = q(sd, .25), sd_q75 = q(sd, .75),
            nb_disp_median = median(nb_disp, na.rm = TRUE),
            cv2_gene_means = var(mean) / mean(mean)^2,
            .groups = 'drop')
diag_summary$nb_disp_median[diag_summary$Scale != 'Raw counts'] <- NA   # not meaningful on log scale
diag_summary$cv2_gene_means[diag_summary$Scale != 'Raw counts'] <- NA
print(as.data.frame(diag_summary), digits = 3)
write.csv(diag_genes,   paste0('Results/DatasetDiagnostics_genes_', tag, '.csv'),   row.names = FALSE)
write.csv(diag_summary, paste0('Results/DatasetDiagnostics_summary_', tag, '.csv'), row.names = FALSE)

# Per-gene distributions (raw counts), one panel per statistic
pDiag <- diag_genes %>%
  filter(Scale == 'Raw counts') %>%
  transmute(Dataset, `Zero fraction` = zero_frac,
            `log10 Mean` = log10(mean),
            `log10 SD` = log10(sd),
            `NB dispersion` = nb_disp) %>%
  tidyr::pivot_longer(-Dataset, names_to = 'Statistic', values_to = 'Value') %>%
  mutate(Statistic = factor(Statistic, levels = c('Zero fraction', 'log10 Mean', 'log10 SD',
                                                  'NB dispersion'))) %>%
  ggplot(aes(Dataset, Value, fill = Dataset)) +
  geom_violin(scale = 'width', linewidth = 0.3, na.rm = TRUE) +
  geom_boxplot(width = 0.15, outliers = FALSE, fill = 'white', linewidth = 0.3, na.rm = TRUE) +
  facet_wrap(~ Statistic, scales = 'free_y', nrow = 1) +
  labs(x = NULL, y = NULL, title = sprintf('Per-gene statistics, raw counts (%d-gene panel)',
                                           length(commonGenes))) +
  theme_classic() + theme(legend.position = 'none', strip.background = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'))
png(paste0('figures/DatasetDiagnostics_', tag, '.png'), height = 3.5, width = 9, units = 'in', res = 300)
print(pDiag); dev.off()

if (use_saved_cutoff && file.exists(cut_file)) {
  cc <- read.csv(cut_file)
  cutoffs_Xi <- setNames(cc$Value, cc$Key)
  message('Loaded saved dependency cutoff: ', cut_file)
} else {
  warning('Saved cutoff not found (', cut_file, '); re-estimating with ', xi_alg)
  dfCount1 <- getCounts(refDataset, n = if (ref_all_cells) NULL else nCells)
  gt1 <- gt[colnames(dfCount1), colnames(dfCount1)]
  cutoffs_Xi <- CalcCut(dfCount1, alpha = alpha, niter = 50, nsamp = nCells, M = xi_M,
                        algorithm = xi_alg, val_mat = gt1, nCores = learn_cores, symmetric = TRUE)
}
# NB: the fallback uses nCores = 4 vs 12 in the cutoff-selection script; with set.seed()
# upstream this should not matter, but prefer the saved file for exact reproducibility.
print(cutoffs_Xi)
DepCut <- as.numeric(cutoffs_Xi[cut_types[[dep_crit]]])
if (!is.finite(DepCut)) stop('Dependency cutoff infeasible for ', dep_crit, '; use dep_crit = "FPR".')

#---------------------
# Directed Xi on GSM8331607: xi[j, k] = xi(j -> k)
#---------------------
dfCount <- getCounts(targetDataset)
gt <- gt[colnames(dfCount), colnames(dfCount)]   # align GT to data gene order
genes <- colnames(dfCount)

# Directed ChatterjeeKNN matrix; no |.|/max-symmetrisation here (that would erase direction).
# Same seed as computeNet() in the cutoff-selection script.
set.seed(0)
xi <- as.matrix(if (xi_alg == 'Xi_KNN') XIM_Mat(dfCount, xi_M, verbose = FALSE) else
  xicor_mpar(dfCount, nCores = learn_cores))
dimnames(xi) <- list(genes, genes)
xi[is.na(xi)] <- 0
diag(xi) <- 0

# Orientation check: XIM_Mat's convention (row->col or col->row) is not documented in the scripts.
# y = x^2 on symmetric x is a function of x but not vice versa, so xi(x->y) >> xi(y->x).
set.seed(1); xt <- rnorm(2000); yt <- xt^2 + rnorm(2000, sd = 0.05)
chk <- as.matrix(if (xi_alg == 'Xi_KNN') XIM_Mat(cbind(x = xt, y = yt), xi_M, verbose = FALSE) else
  xicor_mpar(cbind(x = xt, y = yt), nCores = 1))
if (chk[1, 2] > chk[2, 1]) {
  message('XIM_Mat orientation: m[j,k] = xi(j -> k)  (as assumed)')
} else {
  message('XIM_Mat orientation: m[j,k] = xi(k -> j); transposing so xi[j,k] = xi(j -> k)')
  xi <- t(xi)
}
xiSym <- pmax(xi, t(xi))

#---------------------
# Asymmetric-association cutoff C(n) = c1 * n^-c2
#   'learned': Algorithm 2 (learnCutoffCurve.R) on the TARGET dataset itself at n = nCells.
#              Pairs = target's dependent pairs; N0 = CORUM co-complex pairs (asyNull()).
#              Only null pairs set the cutoff, so the noise level is the target's own. Cached on disk.
#   'formula': fitted on simulations (CutoffSelect_Directed_multiple_KNN.R); uses asy_crit.
#---------------------
if (asy_cut_method == 'formula') {
  fit_file <- paste0('Results/Benchmark/Cutoff_function_fits_FPR_FDR_', tag, '.csv')
  # No silent fallback: the old constants were fitted with the original xi, not KNN.
  if (!file.exists(fit_file)) stop('Run CutoffSelect_Directed_multiple_KNN.R first: ', fit_file)
  fits <- read.csv(fit_file)
  f <- fits[fits$Criterion == asy_crit, ]
  AsyAssCut <- f$c1 * nCells^(-f$c2)
} else {
  source('code/learnCutoffCurve.R')
  curve_file <- paste0('Results/RealData_LearnedCutoff_CORUM_', targetDataset, '_', tag, '.rds')
  curve <- if (file.exists(curve_file)) readRDS(curve_file)
  if (!is.null(curve$cutoff) && identical(curve$settings$n, nCells)) {
    message('Loaded learned cutoff: ', curve_file)
  } else {
    pool        <- getCounts(targetDataset, n = NULL)                  # ALL target cells
    if (nrow(pool) <= nCells) warning('Target has <= nCells cells: Algorithm 2 will bootstrap (ties bias xi)')
    ut_dep      <- which(upper.tri(xiSym) & xiSym > DepCut, arr.ind = TRUE)   # target's dependent pairs
    curve_pairs <- data.frame(g1 = genes[ut_dep[, 1]], g2 = genes[ut_dep[, 2]])
    curve_null  <- asyNull(curve_pairs$g1, curve_pairs$g2)
    message(sprintf('Algorithm 2 on %s: %d cells, %d dependent pairs, %d in N0 (CORUM), n = %d',
                    targetDataset, nrow(pool), nrow(curve_pairs), sum(curve_null), nCells))
    if (sum(curve_null) < 50) warning('Few N0 pairs (', sum(curve_null), '): 95th percentile will be unstable')

    cl <- parallel::makeCluster(learn_cores)
    parallel::clusterEvalQ(cl, {
      library(XICOR)
      source('code/XICOR_mod.R'); source('code/XIMCOR_mod.R')         # XIM_Mat on the workers
    })
    curve <- tryCatch(LearnCutoff(pool, curve_pairs, null = curve_null, n = nCells, B = learn_B,
                                  alpha = alpha, xi_alg = xi_alg, M = xi_M, cl = cl),  # same statistic as XIM_Mat above
                      finally = parallel::stopCluster(cl))
    saveRDS(curve, curve_file)
  }
  AsyAssCut <- curve$cutoff
}
message(sprintf('Dependency cutoff (%s) = %.4f | Asymmetric-association cutoff (%s, n = %d) = %.4f',
                dep_crit, DepCut, if (asy_cut_method == 'learned') 'learned/FPR' else paste0('formula/', asy_crit),
                nCells, AsyAssCut))

#---------------------
# Step 1: dependency test (5) on unique pairs
# Step 2: asymmetric-association test (8) on dependent pairs only
#---------------------
ut  <- which(upper.tri(xi), arr.ind = TRUE)
edges <- data.frame(g1 = genes[ut[, 1]], g2 = genes[ut[, 2]],
                    xi12 = xi[ut], xi21 = t(xi)[ut], xiSym = xiSym[ut],
                    GT = gt[ut]) %>%
  filter(xiSym > DepCut) %>%
  mutate(delta     = xi12 - xi21,
         Asymmetric    = abs(delta) > AsyAssCut,
         Regulator = ifelse(delta >= 0, g1, g2),
         Target    = ifelse(delta >= 0, g2, g1),
         STRING    = case_when(GT %in% TRUE  ~ 'Direct',
                               GT %in% FALSE ~ 'Null',
                               TRUE          ~ 'Indirect/uncertain'),
         TF_fwd    = paste0(Regulator, '_', Target) %in% tf_pairs,   # inferred direction matches TFLink
         TF_rev    = paste0(Target, '_', Regulator) %in% tf_pairs)   # opposite direction in TFLink

message(sprintf('%d dependent pairs; %d Asymmetricly regulated', nrow(edges), sum(edges$Asymmetric)))

# Realised FPR of the asymmetric-association test: share of CORUM null pairs (N0) called asymmetric.
# With asy_cut_method = 'learned' this is in-sample (same N0 set the cutoff), so expect ~alpha.
isNull <- asyNull(edges$g1, edges$g2)
message(sprintf('Asymmetric test on %s: FPR (CORUM N0) = %.4f (%d / %d null pairs), target alpha = %.2f',
                targetDataset, mean(edges$Asymmetric[isNull]), sum(edges$Asymmetric[isNull]),
                sum(isNull), alpha))

# Realised error on the target dataset (non-excluded pairs only)
evalDep <- edges %>% filter(STRING != 'Indirect/uncertain')
nNull   <- sum(gt[upper.tri(gt)] %in% FALSE)
message(sprintf('Dependency test on %s: FP/null = %.4f, FDR (vs STRING) = %.3f, TP = %d',
                targetDataset, sum(evalDep$STRING == 'Null') / nNull,
                mean(evalDep$STRING == 'Null'), sum(evalDep$STRING == 'Direct')))

asymmetric_edges <- edges %>% filter(Asymmetric) %>% arrange(desc(abs(delta)))
print(asymmetric_edges %>% select(Regulator, Target, xi12, xi21, delta, STRING))

# Known TF-target links among asymmetric edges (as in the original script)
tfSub <- asymmetric_edges %>% filter(TF_fwd | TF_rev)
print(tfSub %>% select(Regulator, Target, delta, TF_fwd, TF_rev))
paste0('Correct TF identification = ', round(mean(tfSub$TF_fwd) * 100, 2), '%')

write.csv(edges,        paste0('Results/DependentEdges_',  targetDataset, '_', tag, '.csv'), row.names = FALSE)
write.csv(asymmetric_edges, paste0('Results/AsymmetricAssEdges_',  targetDataset, '_', tag, '.csv'), row.names = FALSE)

#---------------------
# Network of all dependent pairs; Asymmetric associations drawn as arrows
#---------------------
edges <- edges %>%
  mutate(EdgeCategory = factor(ifelse(STRING == 'Direct', 'STRING-validated', 'Dependent only'),
                               levels = c('Dependent only', 'STRING-validated'))) %>%
  arrange(EdgeCategory)                      # draw STRING-validated edges last (on top)

netAll <- graph_from_data_frame(edges %>% select(g1, g2, EdgeCategory), directed = FALSE)
set.seed(0)
layAll <- create_layout(netAll, layout = 'fr')   # node coordinates, reused for the directed GRN
pAll <- ggraph(layAll) +
  geom_edge_link(aes(colour = EdgeCategory), edge_width = 0.4) +
  geom_node_point(size = 1, colour = alpha('magenta4', 0.8), shape = 16) +
  scale_edge_colour_manual(name = 'Edge type', drop = FALSE,
                           values = c('Dependent only' = alpha('grey75', 0.4),
                                      'STRING-validated' = alpha('dodgerblue3', 0.8))) +
  guides(edge_colour = guide_legend(override.aes = list(edge_width = 1), )) +
  theme_graph() +
  theme(legend.position = 'bottom')
png(paste0('figures/Network_', targetDataset, '_AllSignificant_', tag, '.png'),
    height = 8, width = 8, res = 300, units = 'in')
print(pAll)
dev.off()

#---------------------
# Directed GRN (Asymmetric associations only)
#---------------------
if (nrow(asymmetric_edges)) {
  actGenesAll <- unique(c(asymmetric_edges$Regulator, asymmetric_edges$Target))
  netAct <- graph_from_data_frame(asymmetric_edges %>% select(Regulator, Target, STRING),
                                  directed = TRUE,
                                  vertices = data.frame(name = layAll$name))   # all genes, same order as layAll
  pAct <- ggraph(netAct, layout = 'manual', x = layAll$x, y = layAll$y) +
    geom_edge_link(arrow = arrow(length = unit(0.05, 'in'), type = 'closed'),
                   edge_width = 0.4, colour = alpha('firebrick', 0.8),
                   end_cap = circle(0.04, 'in'), start_cap = circle(0.02, 'in')) +
    geom_node_point(size = 1, colour = alpha('magenta4', 0.8), shape = 16) +
    geom_node_text(aes(label = ifelse(name %in% actGenesAll, name, '')),
                   size = 3, repel = TRUE, max.overlaps = 3) +
    theme_graph()
  png(paste0('figures/Network_', targetDataset, '_AsymmetricAss_', tag, '.png'),
      height = 8, width = 8, res = 300, units = 'in')
  print(pAct)
  dev.off()

  #---------------------
  # Ranked-expression scatter for the top asymmetric edge (both orientations)
  #---------------------
  rankDF <- function(g) as.data.frame(apply(dfCount[, g, drop = FALSE], 2, rank,
                                            ties.method = 'random'))
  g1 <- asymmetric_edges$Regulator[1]; g2 <- asymmetric_edges$Target[1]
  set.seed(0); pairDF <- rankDF(c(g1, g2))
  pairDF$celltype <- soPBMC$celltype[rownames(pairDF)]
  sc <- function(x, y, ttl = NULL) ggplot(pairDF, aes(.data[[x]], .data[[y]])) +
    geom_point(aes(color=celltype), size = 1, shape = 16) +
    geom_smooth(method = 'loess', se = TRUE, colour = 'magenta4', fill = alpha('magenta4', 0.2)) +
    theme_classic() + labs(x = paste0(x, ' (Ranked expr.)'), y = paste0(y, ' (Ranked expr.)'),
                           title = ttl)
  pPair <- sc(g1, g2, sprintf('%s → %s\n(ξ = %.3f  vs  ξ reverse = %.3f)',
                              g1, g2, xi[g1, g2], xi[g2, g1])) + sc(g2, g1)
  pPair <- pPair + plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'bottom')
  png(paste0('figures/Scatter_', targetDataset, '_', g1, '_', g2, '_', tag, '.png'),
      height = 4, width = 7, units = 'in', res = 300)
  print(pPair)
  dev.off()

  #---------------------
  # Scatterplot matrix of all genes in asymmetric edges (Supp. Fig. 4), derived not hard-coded
  #---------------------
  actGenes <- unique(c(asymmetric_edges[1:10,]$Regulator, asymmetric_edges[1:10,]$Target))
  set.seed(0); pairsDF <- rankDF(actGenes)
  g <- colnames(pairsDF)

  # Symmetric flag matrix: TRUE where the pair is an Asymmetric association
  act_mat <- matrix(FALSE, length(g), length(g), dimnames = list(g, g))
  for (k in seq_len(nrow(asymmetric_edges))) {
    r <- asymmetric_edges$Regulator[k]; t <- asymmetric_edges$Target[k]
    if (r %in% g && t %in% g) act_mat[r, t] <- act_mat[t, r] <- TRUE
  }

  # par('mfg') gives the (row, col) of the panel being drawn, so no gene names are needed
  hl_panel <- function(x, y, ...) {
    ij <- par('mfg')
    if (act_mat[ij[1], ij[2]]) {
      u <- par('usr')
      rect(u[1], u[3], u[2], u[4], border = 'firebrick', lwd = 5)
    }
    panel.smooth(x, y, pch = 16, cex = 0.2, col = alpha('black', 0.5),
                 col.smooth = 'magenta4', lwd = 1.5)
  }
  pairs(pairsDF, lower.panel = hl_panel, upper.panel = hl_panel)

  png(paste0('figures/Scatter_', targetDataset, '_AsymmetricGenes_', tag, '.png'),
      height = max(6, 1.1 * length(actGenes)), width = max(6, 1.1 * length(actGenes)),
      units = 'in', res = 300)
  pairs(pairsDF, lower.panel = hl_panel, upper.panel = hl_panel)
  dev.off()
}

#---------------------
# TFLink network: asymmetric associations supported by TFLink (either direction).
# Arrow = inferred Regulator -> Target; correct detections (inferred direction matches
# a TFLink TF -> target pair) are marked in green, reversed ones in grey.
#---------------------
tfNet <- asymmetric_edges %>%
  filter(TF_fwd | TF_rev) %>%
  mutate(Detection = factor(ifelse(TF_fwd, 'Correct direction', 'Reversed direction'),
                            levels = c('Reversed direction', 'Correct direction')))  # correct drawn last
message(sprintf('TFLink-supported asymmetric edges: %d | correct direction: %d (%.1f%%)',
                nrow(tfNet), sum(tfNet$TF_fwd), 100 * mean(tfNet$TF_fwd)))

if (nrow(tfNet)) {
  netTF <- graph_from_data_frame(tfNet %>% arrange(Detection) %>%
                                   select(Regulator, Target, Detection), directed = TRUE)
  set.seed(0)
  pTF <- ggraph(netTF, layout = 'fr') +
    geom_edge_link(aes(colour = Detection, width = Detection),
                   arrow = arrow(length = unit(0.08, 'in'), type = 'closed'),
                   end_cap = circle(0.05, 'in'), start_cap = circle(0.03, 'in'),
                   edge_width = 0.6) +
    geom_node_point(size = 2, colour = alpha('magenta4', 0.8), shape = 16) +
    geom_node_text(aes(label = name), size = 3.5, repel = TRUE, max.overlaps = 5) +
    scale_edge_colour_manual(name = 'TFLink', drop = FALSE,
                             values = c('Correct direction'  = alpha('firebrick', 0.5),
                                        'Reversed direction' = alpha('grey50', 0.5))) +
    scale_edge_width_manual(values = c('Correct direction' = 1, 'Reversed direction' = 0.5),
                            guide = 'none') +
    theme_graph()
  png(paste0('figures/Network_', targetDataset, '_TFLink_', tag, '.png'),
      height = 5, width = 6, res = 300, units = 'in')
  print(pTF)
  dev.off()
}
